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
    if x_val is None:
        x_val = np.arange(image.shape[1])
    if y_val is None:
        y_val = np.arange(image.shape[0])

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
