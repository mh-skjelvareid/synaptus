from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
from numpy.typing import NDArray
from scipy.io import loadmat


@dataclass
class UltrasoundDataset:
    raw_data: np.ndarray  # Ultrasound data (e.g., p(t, x))
    fs: float  # Sampling frequency (Hz)
    x_step: float | None = None  # Spatial step along x-axis (meters)
    y_step: float | None = None  # Spatial step along y-axis (meters)
    t_delay: float = 0  # Time delay before recording begins (seconds)
    f_low: float = 0  # Lower frequency limit (Hz)
    f_high: float = float("inf")  # Upper frequency limit (Hz)
    wave_vel: tuple[float] = (1480.0,)  # Wave velocities in each medium (m/s)
    layer_thick: tuple[float] = (float("inf"),)  # Thickness of each layer (meters)


def get_dataset_path(name: str) -> Path:
    """
    Returns the absolute path to a dataset by name.
    Assumes 'datasets/' is located at the same level as the 'python/' directory.
    """
    repo_root = Path(__file__).resolve().parents[2]
    dataset_path = repo_root / "datasets" / name

    if not dataset_path.exists():
        raise FileNotFoundError(f"Dataset not found: {dataset_path}")

    return dataset_path


def load_mat_dataset(filename: str):
    """
    Loads a .mat dataset from the shared datasets folder using scipy.io.
    """
    path = get_dataset_path(filename)
    mat_data = loadmat(path)

    # Load essential fields from the dataset
    if "ptx" in mat_data:
        raw_data = mat_data["ptx"]
    elif "ptxy" in mat_data:
        raw_data = mat_data["ptxy"]
    else:
        raise KeyError("Ultrasound data 'ptx' or 'ptxy' not found in dataset.")

    if "fs" in mat_data:
        fs = mat_data["fs"].item()  # Sampling frequency
    else:
        raise KeyError("Sampling frequency 'fs' not found in dataset.")

    # Create ultrasound dataset instance
    dataset = UltrasoundDataset(raw_data=raw_data, fs=fs)

    # Load optional fields
    if "xStep" in mat_data:
        dataset.x_step = mat_data["xStep"].item()  # Spatial step size (x-axis)
    if "yStep" in mat_data:
        dataset.y_step = mat_data["yStep"].item()
    if "tDelay" in mat_data:
        dataset.t_delay = mat_data["tDelay"].item()
    if "cc" in mat_data:
        dataset.wave_vel = tuple(np.ravel(mat_data["cc"]))  # Wave velocities in each medium
    if "thick" in mat_data:
        dataset.layer_thick = tuple(np.ravel(mat_data["thick"]))
    if "fLow" in mat_data:
        dataset.f_low = mat_data["fLow"].item()
    if "fHigh" in mat_data:
        dataset.f_high = mat_data["fHigh"].item()

    return dataset
