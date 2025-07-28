from math import ceil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.colors import Colormap
from numpy.typing import NDArray
from scipy.signal import hilbert


def nextpow2(n: float) -> int:
    """Round up number to next power of 2"""
    return int(2 ** ceil(np.log2(n)))


def get_envelope(arr: NDArray, axis: int = 0) -> NDArray:
    """Calculate envelope of signal using Hilbert transform

    Parameters
    ----------
    arr : NDArray
        Vector or tensor containing signal(s) in time-domain.
    axis : int, optional
        Axis corresponding to time dimension of arr, by default 0

    Returns
    -------
    NDArray
        Array with envelope values, same shape as arr.
    """
    return np.abs(hilbert(arr, axis=axis))  # type: ignore


def log_image(im: NDArray) -> NDArray:
    """Normalize image and transform to deciBel (log) scale"""
    return 20 * np.log10(im / np.max(im))


def plot_us_image(
    image: NDArray,
    axes: Axes | None = None,
    x_val: NDArray | None = None,
    y_val: NDArray | None = None,
    x_label: str = "",
    y_label: str = "",
    title: str = "",
    min_db: float = -60,
    max_db: float = 0,
    figsize: tuple[float, float] = (6, 4),
    cmap: str | Colormap = "viridis",
) -> None:
    """Plot ultrasound image (abs.val.) on logarithmic scale, with colorbar

    Parameters
    ----------
    image : NDArray
        Ultrasound image, typically array of float or complex
    axes : Axes | None, optional
        Axes into which to plot. If None, a new figure and Axes object is created.
    x_val : NDArray | None, optional
        Vector of values corresponding to horizontal axis (columns) of image.
    y_val : NDArray | None, optional
        Vector of values corresponding to vertical axis (rows) of image.
    x_label : str, optional
        Label displayed below X axis in plot, by default ""
    y_label : str, optional
        Label displayed beside Y axis in plot, by default ""
    title : str, optional
        Title displayed above plot, by default ""
    min_db : float, optional
        Minimum value displayed (in dB), by default -60
    max_db : float, optional
        Maximum value displayed (in dB), by default 0
        Since the image is normalized so that the maximum value corresponds to 0 dB,
        there is usually no need to change this.
    figsize : tuple[float, float], optional
        Size of figure (if axes parameter is None), by default (6, 4)
    cmap: str | Colormap
        Colormap used when plotting image, by default "viridis"

    Raises
    ------
    ValueError
        Raises error if image array is not 2D.
    """
    if image.ndim != 2:
        raise ValueError("Image must be a 2D array")
    if axes is None:
        _, axes = plt.subplots(figsize=figsize)
    x_val = x_val if x_val is not None else np.arange(image.shape[1])
    y_val = y_val if y_val is not None else np.arange(image.shape[0])

    # Plot image
    im_handle = axes.imshow(
        log_image(np.abs(image)),
        interpolation="none",
        extent=(
            x_val[0],
            x_val[-1],
            y_val[-1],
            y_val[0],
        ),
        aspect="auto",
        vmin=min_db,
        vmax=max_db,
        cmap=cmap,
    )

    # Create colorbar
    cbar = plt.colorbar(im_handle)

    # Set text labels
    axes.set_xlabel(x_label)
    axes.set_ylabel(y_label)
    axes.set_title(title)
    cbar.ax.set_ylabel("dB", rotation=270, labelpad=15)

    # Show plot
    plt.show()


def make_omega_vec(n_fft: int, fs: float) -> NDArray:
    """Make vector of omega values corresponding to FFT bins.

    Parameters
    ----------
    n_fft : int
        Number of points in FFT
    fs : float
        Sampling frequency in Hz

    Returns
    -------
    NDArray
        Vector of omega (angular frequency) values.
    """
    return np.fft.fftshift(np.fft.fftfreq(n_fft, 1 / fs)) * (2 * np.pi)


def make_k_vec(n_fft: int, step: float) -> NDArray:
    """Make vector of k values corresponding to FFT bins.

    Parameters
    ----------
    n_fft : int
        Number of points in FFT
    step : float
        Spatial step size (in meters) corresponding to FFT bins

    Returns
    -------
    NDArray
        Vector of k (spatial frequency) values.
    """
    return np.fft.fftshift(np.fft.fftfreq(n_fft, step)) * (2 * np.pi)


def calc_omega_passband(omega: NDArray, f_low, f_high) -> NDArray:
    """Calculate passband of omega values based on frequency limits.

    Parameters
    ----------
    omega : NDArray
        Vector of omega values.
    f_low : float
        Lower frequency limit in Hz.
    f_high : float
        Upper frequency limit in Hz.

    Returns
    -------
    NDArray
        Vector of boolean values indicating passband.

    """
    return (omega >= (2 * np.pi * f_low)) & (omega <= (2 * np.pi * f_high))


def calc_depth_resolution(wave_velocity: float, f_low: float, f_high: float) -> float:
    """Calculate depth resolution based on wave velocity and transducer bandwidth.

    Parameters
    ----------
    wave_velocity : float
        Wave velocity in the medium (in m/s).
    f_low : float
        Lower frequency limit of the transducer (in Hz).
    f_high : float
        Upper frequency limit of the transducer (in Hz).

    Returns
    -------
    float
        Depth resolution (in meters).
    """
    return (wave_velocity / 2) / (f_high - f_low)


def make_coord_grids(*args) -> tuple[NDArray, ...]:
    """Make coordinate grids based on coordinate vectors

    Wrapper for numpy.meshgrid using `indexing='ij'` to ensure
    that the first dimension corresponds to the first coordinate vector.

    Parameters
    ----------
    *args : tuple
        Coordinate vectors for each dimension, e.g. (omega_vec, kx_vec) for 2D data

    Returns
    -------
    list[NDArray]
        Tuple containing grids for each dimension, e.g. 2D matrices omega_mat and kx_mat
        for the 2D case of (omega,kx) coordinates.

    See also:
    --------
    numpy.meshgrid
    """
    return np.meshgrid(*args, indexing="ij")


def make_kz_grid(
    wave_velocity: float, omega_grid: NDArray, *k_grids: NDArray
) -> tuple[NDArray, NDArray]:
    """Calculate wavenumber kz for every combination of omega and kx

    Parameters
    ----------
    wave_velocity : float
        Wave velocity of medium in m/s.
    omega_grid : NDArray
        _description_
    k_grids : list[NDArray]
        List of 1 or 2 arrays containing kx and ky values.
        If only one array is provided, it is assumed to be kx.
        If two arrays are provided, they are assumed to be kx and ky.

    Returns
    -------
    tuple[NDArray, NDArray]
        Grid of kz values and boolean array indicating real wave index.

    Raises
    ------
    ValueError
        If k_grids contains more than 2 elements, or if it is empty.
    """

    n_k_grids = len(k_grids)
    if n_k_grids == 1:
        kx_grid = k_grids[0]
        kz_grid_sq = ((2 / wave_velocity) ** 2) * (omega_grid**2) - kx_grid**2
    elif n_k_grids == 2:
        kx_grid, ky_grid = k_grids
        kz_grid_sq = ((2 / wave_velocity) ** 2) * (omega_grid**2) - kx_grid**2 - ky_grid**2
    else:
        raise ValueError("k_grids must contain 1 element (kx) or 2 elements (kx, ky)")

    real_wave_index = kz_grid_sq >= 0
    KZ = np.sqrt(kz_grid_sq * real_wave_index)
    return KZ, real_wave_index
