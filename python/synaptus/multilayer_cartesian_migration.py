from pathlib import Path

import numpy as np
from docstring_inheritance import NumpyDocstringInheritanceMeta
from numpy.typing import NDArray
from rich import print
from scipy.io import loadmat
from utils import (
    calc_depth_resolution,
    calc_omega_passband,
    make_coord_grids,
    make_k_vec,
    make_kz_grid,
    make_omega_vec,
    nextpow2,
)


class MultilayerCartesianPulseEchoData(metaclass=NumpyDocstringInheritanceMeta):
    """Class for performing phase shift migration on pulse-echo data."""

    def __init__(
        self,
        raw_data: NDArray,
        fs: float,
        f_low: float,
        f_high: float,
        x_step: float,
        y_step: float | None = None,
        t_delay: float = 0.0,
        wave_velocities: tuple[float] = (1500,),
        layer_thicknesses: tuple[float] | None = None,
        nfft_t: int | None = None,
        nfft_x: int | None = None,
        nfft_y: int | None = None,
    ) -> None:
        """_summary_

        Parameters
        ----------
        raw_data : NDArray
            2D or 3D array containing pulse-echo data,
            The first dimension corresponds to time, the second to spatial dimension x,
            and the third to spatial dimension y (if present).
        fs : float
            Sampling frequency of the data (in time domain) in Hz.
        f_low : float
            Low frequency cutoff for the transducer band in Hz.
        f_high : float
            High frequency cutoff for the transducer band in Hz.
        x_step : float
            Spatial step size in the x direction in meters.
            This is the distance between adjacent transducer positions.
        y_step : float | None, optional
            Spatial step size in the y direction in meters (for 3D data).
        t_delay : float, optional
            Time delay from pulse transmission to start of data acquisition, in seconds.
            By default 0.0
        wave_velocities : tuple[float], optional
            Wave velocities in the medium(s) through which the wave travels, in
            meters per second. If only one value is provided, it is assumed that the
            wave velocity is constant throughout the medium. If not specified, a
            default value of 1500 m/s is used, which is typical for water/soft tissue.
        layer_thicknesses : tuple[float] | None, optional
            Thicknesses of layers through which the wave travels, in meters. If not
            specified, it is assumed that there is only one layer, and the thickness is
            estimated the raw data. If multiple layers are specified, the number of
            layer thicknesses must match the number of wave velocities.
        nfft_t : int | None, optional
            Number of points for FFT in time dimension. If None, set to the next power of 2
            greater than the number of time samples in the raw data.
        nfft_x : int | None, optional
            Number of points for FFT in x dimension. If None, set to the next power of 2
            greater than the number of spatial samples in the x dimension.
        nfft_y : int | None, optional
            Number of points for FFT in y dimension. If None, set to the next power of 2
            greater than the number of spatial samples in the y dimension.

        Notes
        -----
        - When processing data in the Fourier domain, "aliasing" artefacts can appear in
          the focused image. This problem can be mitigated by increasing the number of
          points in the FFTs used, i.e. `nfft_t`, `nfft_x`, or `nfft_y`. In the
          time/spatial domain, this corresponds to zero-padding the data.
        """
        # Input data
        self.raw_data = raw_data
        self.ndim = raw_data.ndim
        if self.ndim not in (2, 3):
            raise ValueError("raw_data must be a 2D or 3D array.")
        self.nt = raw_data.shape[0]
        self.nx = raw_data.shape[1]
        self.ny = raw_data.shape[2] if self.ndim == 3 else None

        # Acquisition parameters
        self.fs = fs
        self.f_low = f_low
        self.f_high = f_high
        self.t_delay = t_delay
        self.x_step = x_step
        self.y_step = y_step
        if self.ndim == 3 and y_step is None:
            raise ValueError("y_step must be specified for 3D data.")

        # Make time-space coordinate vectors
        self.time_vec = np.arange(self.nt) / self.fs + self.t_delay
        self.x_vec = np.arange(self.nx) * self.x_step
        self.y_vec = np.arange(self.ny) * self.y_step if self.ndim == 3 else np.empty(0)  # type:ignore

        # Medium properties
        self.wave_velocities = wave_velocities
        self.layer_thicknesses = layer_thicknesses
        if (len(wave_velocities) > 1) and (len(wave_velocities)) != len(layer_thicknesses):
            raise ValueError(
                "If multiple wave velocities are provided, "
                "the number of layer thicknesses must match."
            )
        if layer_thicknesses is None:
            layer_thicknesses = (
                self.time_vec[-1] * (self.wave_velocities[0] / 2),
            )  # Layer thickness for single layer = end of measurement

        ## FFT settings
        self.nfft_t = nfft_t if nfft_t is not None else nextpow2(self.nt)
        self.nfft_x = nfft_x if nfft_x is not None else nextpow2(self.nx)
        if self.ndim == 3:
            self.nfft_y = nfft_y if nfft_y is not None else nextpow2(self.ny)  # type:ignore
        else:
            self.nfft_y = -1  # Placeholder value

        # Make frequency domain coordinate vectors
        self.omega_vec_full = make_omega_vec(self.nfft_t, self.fs)
        self.omega_passband = calc_omega_passband(self.omega_vec_full, self.f_low, self.f_high)
        self.omega_vec = self.omega_vec_full[self.omega_passband]
        self.omega_vec_full = make_omega_vec(self.nfft_t, self.fs)
        self.kx_vec = make_k_vec(self.nfft_x, self.x_step)
        self.ky_vec = make_k_vec(self.nfft_y, self.y_step) if self.ndim == 3 else np.empty(0)  # type:ignore

        # Make frequency-domain coordinate grids
        if self.ndim == 2:
            self.omega_grid, self.kx_grid = make_coord_grids(self.x_vec, self.y_vec)
            self.ky_grid = np.empty(0)  # Placeholder
        else:
            self.omega_grid, self.kx_grid, self.ky_grid = make_coord_grids(self.x_vec, self.y_vec)

        # Perform Fourier transform on the raw data, and phase shift to t=0
        self.wavefield = self._time_shift_to_t_zero(self._fourier_transform())

    def _fourier_transform(self) -> NDArray:
        """Perform Fourier transform on the raw data."""
        if self.ndim == 2:
            wavefield = np.fft.fftshift(np.fft.fftn(self.raw_data, s=(self.nfft_t, self.nfft_x)))
        else:
            wavefield = np.fft.fftshift(
                np.fft.fftn(self.raw_data, s=(self.nfft_t, self.nfft_x, self.nfft_y))
            )
        return wavefield[self.omega_passband]  # Crop to pos. omega in transducer passband

    def _inverse_fourier_transform(self, wavefield: NDArray) -> NDArray:
        """Perform inverse Fourier transform on a frequency-domain wavefield."""
        if self.ndim == 2:  # 2D
            return np.fft.ifftn(np.fft.ifftshift(wavefield, axes=(1,)))
        else:  # 3D
            return np.fft.ifftn(np.fft.ifftshift(wavefield, axes=(1, 2)))

    def _time_shift_to_t_zero(self, wavefield: NDArray) -> NDArray:
        """Apply negative time shift (as phase shift) to align the wavefield to t=0."""
        return wavefield * np.exp(-1j * self.omega_grid * self.t_delay)

    def z_shift_wavefield(self, wavefield: NDArray, wave_velocity: float, dz: float):
        """Apply a phase shift to the wavefield to account for a depth shift.

        Parameters
        ----------
        wavefield : NDArray
            Wavefield in the frequency domain.
        wave_velocity : float
            Wave velocity in the medium (in m/s).
        dz : float
            Depth shift to apply (in meters).

        Returns
        -------
        NDArray
            Wavefield with the applied depth shift.
        """
        # Calculate the kz grid based on the wave velocity and the frequency-domain grids
        if self.ndim == 2:
            kz_grid, real_wave_index = make_kz_grid(wave_velocity, self.omega_grid, self.kx_grid)
        else:
            kz_grid, real_wave_index = make_kz_grid(
                wave_velocity, self.omega_grid, self.kx_grid, self.ky_grid
            )
        return wavefield * np.exp(1j * kz_grid * dz) * real_wave_index


class PhaseShiftMigration(MultilayerCartesianPulseEchoData):
    """Class for performing phase shift migration on pulse-echo data."""

    def __init__(self, *args, **kwargs) -> None:
        """Initialize the PhaseShiftMigration class."""
        super().__init__(*args, **kwargs)

    def calc_phase_shift_tensor(self, wave_velocity: float) -> tuple[NDArray, NDArray]:
        # Make kz grid
        if self.ndim == 2:
            kz_grid, real_wave_index = make_kz_grid(wave_velocity, self.omega_grid, self.kx_grid)
        else:
            kz_grid, real_wave_index = make_kz_grid(
                wave_velocity, self.omega_grid, self.kx_grid, self.ky_grid
            )

        # Calculate depth resolution
        dz = calc_depth_resolution(wave_velocity, self.f_low, self.f_high)

        return np.exp(1j * kz_grid * dz), real_wave_index

    def phase_shift_migrate(self) -> list[NDArray]:
        """Perform phase shift migration on the wavefield."""
        wavefield = self.wavefield.copy()
        images = []

        for wave_velocity, layer_thickness in zip(self.wave_velocities, self.layer_thicknesses):
            # Get reolution and number of depth samples in current layer
            dz = calc_depth_resolution(wave_velocity, self.f_low, self.f_high)
            n_depth_samples = int(layer_thickness / dz)
            layer_z_vec = np.arange(n_depth_samples) * dz

            # Preallocate array for fucused image in this layer
            layer_image = np.zeros(
                shape=(n_depth_samples,) + wavefield.shape[1:], dtype=np.complex128
            )

            # Calculate phase shift tensor for this wave velocity
            phase_shift_tensor, real_wave_index = self.calc_phase_shift_tensor(wave_velocity)
            wavefield *= real_wave_index  # Apply real wave index to remove non-physical components

            # Phase shift line by line
            for line_ind in range(n_depth_samples):
                layer_image[line_ind] = np.fft.ifft(
                    np.sum(wavefield, axis=0, keepdims=True), axis=1
                )
                wavefield *= phase_shift_tensor

            # TODO: Fix potential round-off error at layer interface(?)

            # Save image for this layer (absolute value, without zero-padding)
            if self.ndim == 2:
                layer_image = np.abs(layer_image[:, : self.nx])
            else:
                layer_image = np.abs(layer_image[:, : self.nx, : self.ny])
            images.append(layer_image)

        return images


if __name__ == "__main__":
    example_data_path = Path().resolve().parent.parent / "datasets" / "LineScan2D_WireTargets.mat"
    example_data = loadmat(example_data_path)
    raw_data = example_data["ptx"]  # 2D ultrasound data p(t,x)
    fs = example_data["fs"].item()  # Sampling frequency
    x_step = example_data["xStep"].item()  # Spatial step size (x-axis)
    sound_vel = example_data["cc"].item()  # Sound velocity
    t_delay = example_data["tDelay"].item()  # Pulse recording delay

    f_low = 0.4e6  # Lower cutoff freq., transducer band
    f_high = 2.5e6  # Upper cutoff freq., transducer band

    dataset = MultilayerCartesianPulseEchoData(
        raw_data=raw_data,
        fs=fs,
        f_low=f_low,
        f_high=f_high,
        x_step=x_step,
        t_delay=t_delay,
        wave_velocities=(sound_vel,),
    )
    print(vars(dataset))
