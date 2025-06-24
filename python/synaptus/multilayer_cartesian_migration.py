from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from utils import (
    calc_depth_resolution,
    calc_omega_passband,
    make_coord_grids,
    make_k_vec,
    make_omega_vec,
    nextpow2,
)


class MultilayerCartesianPulseEchoData:
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
        sound_velocities: tuple[float] = (1500,),
        layer_thicknesses: tuple[float] | None = None,
        nfft_t: int | None = None,
        nfft_x: int | None = None,
        nfft_y: int | None = None,
    ) -> None:
        """_summary_

        Parameters
        ----------
        raw_data : NDArray
            2D or 3D array containing ultrasound data,
            The first dimension corresponds to time, the second to spatial dimension x,
            and the third to spatial dimension y (if present).
        fs : float
            Sampling frequency of the ultrasound data in Hz.
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
        sound_velocities : tuple[float], optional
            Sound velocities in the medium(s) through which the ultrasound travels, in
            meters per second. If only one value is provided, it is assumed that the
            sound velocity is constant throughout the medium. If not specified, a
            default value of 1500 m/s is used, which is typical for water/soft tissue.
        layer_thicknesses : tuple[float] | None, optional
            Thicknesses of layers through which the ultrasound
            travels, in meters. If not specified, it is assumed that there is only one
            layer, and the thickness is estimated the raw data. If multiple layers are
            specified, the number of layer thicknesses must match the number of
            sound velocities.
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

        # Medium properties
        self.sound_velocities = sound_velocities
        self.layer_thicknesses = layer_thicknesses
        if (len(sound_velocities) > 1) and (len(sound_velocities)) != len(layer_thicknesses):
            raise ValueError(
                "If multiple sound velocities are provided, "
                "the number of layer thicknesses must match."
            )
        if layer_thicknesses is None:
            layer_thicknesses = (
                self.time_vec[-1] * (self.sound_velocities[0] / 2),
            )  # Layer thickness for single layer = end of measurement

        ## FFT settings
        self.nfft_t = nfft_t if nfft_t is not None else nextpow2(self.nt)
        self.nfft_x = nfft_x if nfft_x is not None else nextpow2(self.nx)
        if self.ndim == 3:
            self.nfft_y = nfft_y if nfft_y is not None else nextpow2(self.ny)  # type:ignore
        else:
            self.nfft_y = None

        # Make time-space coordinate vectors
        self.time_vec = np.arange(self.nt) / self.fs + self.t_delay
        self.x_vec = np.arange(self.nx) * self.x_step
        self.y_vec = np.arange(self.ny) * self.y_step if self.ndim == 3 else np.empty(0)  # type:ignore

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
