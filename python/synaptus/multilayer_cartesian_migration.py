from pathlib import Path

import numpy as np
from docstring_inheritance import NumpyDocstringInheritanceMeta
from numpy.typing import NDArray
from scipy.interpolate import RegularGridInterpolator

from .utils import (
    calc_depth_resolution,
    calc_omega_passband,
    make_coord_grids,
    make_k_vec,
    make_kz_grid,
    make_omega_vec,
    nextpow2,
)


class MultilayerCartesianPulseEchoData(metaclass=NumpyDocstringInheritanceMeta):
    """Dataset class for pulse-echo data in cartesian coordinate system"""

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
        omega_upsampling_factor: int = 1,
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
        omega_upsampling_factor : int, optional
            Factor by which to upsample the omega vector. This is useful for improving
            the resolution of the frequency-domain representation. Default is 1 (no
            upsampling). Ignored if nfft_t is specified.


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
        if (len(wave_velocities) > 1) and (len(wave_velocities)) != len(layer_thicknesses):
            raise ValueError(
                "If multiple wave velocities are provided, "
                "the number of layer thicknesses must match."
            )
        if layer_thicknesses is None:
            layer_thicknesses = (
                self.time_vec[-1] * (self.wave_velocities[0] / 2),
            )  # Layer thickness for single layer = end of measurement
        self.layer_thicknesses = layer_thicknesses

        ## FFT settings
        self.nfft_t = nfft_t if nfft_t is not None else nextpow2(self.nt) * omega_upsampling_factor
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
            self.omega_grid, self.kx_grid = make_coord_grids(self.omega_vec, self.kx_vec)
            self.ky_grid = np.empty(0)  # Placeholder
        else:
            self.omega_grid, self.kx_grid, self.ky_grid = make_coord_grids(
                self.omega_vec, self.kx_vec, self.ky_vec
            )

        # Perform Fourier transform on the raw data, and phase shift to t=0
        self.wavefield = self._time_shift_to_t_zero(self._fourier_transform())

    def _fourier_transform(self) -> NDArray:
        """
        Performs an n-dimensional Fourier transform on the raw input data and returns
        the frequency-domain wavefield cropped to the positive frequency components
        within the transducer passband.

        Returns:
            NDArray: The Fourier-transformed wavefield, restricted to the positive
            frequency passband.

        Notes:
            - For 2D data, the transform is performed over time and x axes (time, x).
            - For 3D data, the transform is performed over time, x, and y axes (time, x,
              y).
            - Only positive (omega) frequencies within the transducer's passband are
              returned.
        """

        if self.ndim == 2:
            wavefield = np.fft.fftshift(
                np.fft.fftn(self.raw_data, s=(self.nfft_t, self.nfft_x), axes=(0, 1))
            )
        else:
            wavefield = np.fft.fftshift(
                np.fft.fftn(
                    self.raw_data, s=(self.nfft_t, self.nfft_x, self.nfft_y), axes=(0, 1, 2)
                )
            )
        return wavefield[self.omega_passband]  # Crop to pos. omega in transducer passband

    def _time_shift_to_t_zero(self, wavefield: NDArray) -> NDArray:
        """
        Apply a negative time shift to the input wavefield to align it with t=0.

        This method multiplies the input wavefield by a complex exponential factor,
        effectively applying a phase shift in the frequency domain. The phase shift
        corresponds to a negative time delay (`-self.t_delay`), which aligns the wavefield
        so that its reference time is zero.

        Parameters
        ----------
        wavefield : NDArray
            The input wavefield array, typically in the frequency domain.

        Returns
        -------
        NDArray
            The time-shifted wavefield, aligned to t=0.

        """
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
        """
        Computes the phase shift tensor for wavefield extrapolation at a specified wave velocity.

        This method constructs the vertical wavenumber (kz) grid based on the spatial and frequency grids,
        then calculates the phase shift tensor used for propagating wavefields in depth. The depth increment (dz)
        is determined from the frequency range and wave velocity. The method supports both 2D and 3D grids.

        Parameters
        ----------
        wave_velocity : float
            The propagation velocity of the wave (in meters per second).

        Returns
        -------
        phase_shift_tensor : NDArray
            Complex-valued tensor representing the phase shift for each grid point.
        real_wave_index : NDArray
            Boolean or integer array indicating the indices of physically meaningful (real) wavenumbers.

        Notes
        -----
        The phase shift tensor is computed as exp(1j * kz * dz), where kz is the vertical wavenumber grid
        and dz is the depth increment. This is commonly used in wavefield extrapolation methods such as
        phase-shift migration.
        """
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

    def phase_shift_migrate(self) -> tuple[list[NDArray], list[np.ndarray]]:
        """
        Perform phase shift migration on the wavefield for each subsurface layer.

        This method applies the phase shift migration algorithm to the current wavefield,
        propagating it through each defined subsurface layer using the corresponding wave velocity
        and layer thickness. For each layer, the wavefield is phase-shifted in the frequency domain,
        and an image is constructed by inverse Fourier transforming the summed wavefield at each depth step.

        Returns:
            list[NDArray]: A list of 2D or 3D NumPy arrays (depending on self.ndim), where each array
                represents the migrated image for a corresponding subsurface layer. The images are
                absolute-valued and cropped to the original spatial dimensions.

        Notes:
            - The method assumes that self.wavefield, self.wave_velocities, self.layer_thicknesses,
              self.f_low, self.f_high, self.nx, self.ny, and self.ndim are properly initialized.
            - Non-physical wave components are suppressed using the real wave index.
        """

        # Copy the wavefield to avoid modifying the original data
        interface_wavefield = self.wavefield.copy()

        # Preallocate list for images (one per layer)
        images = []
        z_vecs = []

        for wave_velocity, layer_thickness in zip(self.wave_velocities, self.layer_thicknesses):
            # Create a copy of the interface wavefield (migrated to top of layer) for this layer
            wavefield = interface_wavefield.copy()

            # Get resolution and number of depth samples in current layer
            dz = calc_depth_resolution(wave_velocity, self.f_low, self.f_high)
            layer_nz = int(layer_thickness / dz)

            # Preallocate array for focused image in this layer
            layer_image = np.zeros(shape=(layer_nz,) + wavefield.shape[1:], dtype=np.complex128)

            # Calculate phase shift tensor for this wave velocity
            phase_shift_tensor, real_wave_index = self.calc_phase_shift_tensor(wave_velocity)
            wavefield *= real_wave_index  # Apply real wave index to remove non-physical components

            # Phase shift line by line
            for line_ind in range(layer_nz):
                layer_image[line_ind] = np.fft.ifft(
                    np.sum(wavefield, axis=0, keepdims=True), axis=1
                )
                wavefield *= phase_shift_tensor

            # Save image for this layer (absolute value, without zero-padding)
            if self.ndim == 2:
                layer_image = np.abs(layer_image[:, : self.nx])
            else:
                layer_image = np.abs(layer_image[:, : self.nx, : self.ny])
            images.append(layer_image)

            # Create and save z coordinate vector for this layer
            layer_z_vec = np.arange(layer_nz) * dz
            z_vecs.append(layer_z_vec)

            # Migrate wavefield to next layer
            interface_wavefield = self.z_shift_wavefield(
                interface_wavefield, wave_velocity, layer_thickness
            )

        return images, z_vecs


class MultilayerOmegaKMigration(MultilayerCartesianPulseEchoData):
    def __init__(self, *args, omega_upsampling_factor=4, **kwargs) -> None:
        """Initialize the StoltMigration class."""
        kwargs["omega_upsampling_factor"] = omega_upsampling_factor
        super().__init__(*args, **kwargs)

    def _calc_layer_kz_vector(self, layer_index) -> NDArray:
        """Calculate the kz vector for a given layer index."""
        relative_bandwidth = (self.f_high - self.f_low) / (self.fs / 2)
        nfft_z = nextpow2(self.nt * relative_bandwidth)
        dz = calc_depth_resolution(self.wave_velocities[layer_index], self.f_low, self.f_high)
        kz_vec = (2 * np.pi) * (
            2 * self.f_low / self.wave_velocities[layer_index] + np.arange(nfft_z) / (dz * nfft_z)
        )
        return kz_vec

    def _stolt_transform(self, wavefield: NDArray, layer_index: int) -> NDArray:
        """
        Applies the Stolt migration transform to the input wavefield for a specified
        layer.

        The Stolt transform is a frequency-wavenumber domain migration technique used in
        seismic imaging. This method interpolates the input wavefield from (omega, kx[,
        ky]) coordinates to new (omega, kx[, ky]) coordinates according to the Stolt
        mapping for the given layer's velocity. It also applies the appropriate
        amplitude scaling factor.

        Parameters:
            wavefield (NDArray): The input wavefield in the frequency-wavenumber domain.
                Shape should match (omega, kx) for 2D or (omega, kx, ky) for 3D.
            layer_index (int): Index of the layer for which to perform the Stolt
            transform.
                Used to select the appropriate velocity and kz vector.

        Returns:
            NDArray: The Stolt-migrated wavefield, interpolated and amplitude-corrected,
                with the same shape as the input wavefield.

        Notes:
            - For 2D data, the transform is performed over (omega, kx).
            - For 3D data, the transform is performed over (omega, kx, ky).
            - Uses linear interpolation and zero-filling for out-of-bounds values.
            - Requires precomputed frequency (omega_vec), wavenumber (kx_vec, [ky_vec]),
              and velocity arrays.
        """

        # Calculate the kz vector for the current layer
        kz_vec = self._calc_layer_kz_vector(layer_index)

        # Create meshgrid coordinate matrices
        if self.ndim == 2:
            KZ_interp, KX_interp = make_coord_grids(kz_vec, self.kx_vec)
            OMEGA_interp = np.sqrt(KZ_interp**2 + KX_interp**2) * (
                self.wave_velocities[layer_index] / 2
            )

            # Calculate amplitude scaling factor
            Akzkx = 1 / (1 + (KX_interp**2) / (KZ_interp**2))

            # Create interpolator object
            interpolator = RegularGridInterpolator(
                points=(self.omega_vec, self.kx_vec),
                values=wavefield,
                bounds_error=False,
                fill_value=0,
                method="linear",
            )

            # Interpolate the wavefield to new omega coordinates
            return interpolator((OMEGA_interp, KX_interp)) * Akzkx

        else:  # ndim == 3
            KZ_interp, KX_interp, KY_interp = make_coord_grids(kz_vec, self.kx_vec, self.ky_vec)
            OMEGA_interp = np.sqrt(KZ_interp**2 + KX_interp**2 + KY_interp**2) * (
                self.wave_velocities[layer_index] / 2
            )

            # Calculate amplitude scaling factor
            Akzkx = 1 / (1 + (KX_interp**2 + KY_interp**2) / (KZ_interp**2))

            # Create interpolator object
            interpolator = RegularGridInterpolator(
                points=(self.omega_vec, self.kx_vec, self.ky_vec),
                values=wavefield,
                bounds_error=False,
                fill_value=0,
                method="linear",
            )

            # Interpolate the wavefield to new omega coordinates
            return interpolator((OMEGA_interp, KX_interp, KY_interp)) * Akzkx

    def mulok_migrate(self) -> tuple[list[NDArray], list[NDArray]]:
        """
        Perform multilayer omega-k (ω-k) migration on the wavefield.

        This method applies Stolt migration sequentially to each subsurface layer,
        as defined by the provided wave velocities and layer thicknesses. For each layer:
          - The depth sampling interval (`dz`) is computed based on the local wave velocity
            and frequency bounds.
          - The wavefield is migrated using a Stolt transform, and the resulting image for
            the current layer is extracted.
          - The wavefield is then propagated (shifted) to the next interface for further migration.

        Returns:
            tuple[list[NDArray], list[NDArray]]: A tuple containing a list of migrated
            images (one per layer) and a list of corresponding z coordinate vectors for
            each layer.

        Notes:
            - The migration is performed in the frequency-wavenumber (ω-k) domain.
        """
        interface_wavefield = self.wavefield.copy()
        images = []
        z_vecs = []

        # Iterate over each layer, creating focused images and migrating between layer interfaces
        for wave_velocity, layer_index in zip(
            self.wave_velocities, range(len(self.layer_thicknesses))
        ):
            # Get resolution and number of depth samples in current layer
            dz = calc_depth_resolution(wave_velocity, self.f_low, self.f_high)
            layer_nz = int(self.layer_thicknesses[layer_index] / dz)

            # Focus image for current layer
            stolt_migrated_wavefield = self._stolt_transform(interface_wavefield, layer_index)
            if self.ndim == 2:
                layer_image = np.fft.ifftn(np.fft.ifftshift(stolt_migrated_wavefield, axes=(1,)))
                layer_image = np.abs(layer_image[:layer_nz, : self.nx])
            else:
                layer_image = np.fft.ifftn(np.fft.ifftshift(stolt_migrated_wavefield, axes=(1, 2)))
                layer_image = np.abs(layer_image[:layer_nz, : self.nx, : self.ny])
            images.append(layer_image)

            # Create and save z coordinate vector for this layer
            layer_z_vec = np.arange(layer_nz) * dz
            z_vecs.append(layer_z_vec)

            # Migrate wavefield to next layer
            interface_wavefield = self.z_shift_wavefield(
                interface_wavefield, wave_velocity, self.layer_thicknesses[layer_index]
            )

        return images, z_vecs


if __name__ == "__main__":
    pass
