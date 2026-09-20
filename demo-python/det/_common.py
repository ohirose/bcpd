"""Shared execution utilities for the DET demos."""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

DEMO = Path(__file__).resolve().parent.parent
ROOT = DEMO.parent
VIEWERS = DEMO / "viewers"
sys.path.insert(0, str(VIEWERS))
from _viewer import compare


def data_file(name):
    for folder in (ROOT / "data", ROOT / "data/local"):
        filename = folder / name
        if filename.exists():
            return filename
    raise FileNotFoundError(f"Data file not found: {name}")


def run(target, source, target_values, source_values, options, trajectory=False):
    """Run DET and return the transformed source domain."""
    bcpd = ROOT / ("win/bcpd.exe" if os.name == "nt" else "bcpd")
    if not bcpd.exists():
        raise FileNotFoundError("BCPD is not built; run make in the repository root")

    with tempfile.TemporaryDirectory(prefix="det-demo-") as work:
        work = Path(work)
        arrays = target, source, target_values, source_values
        names = "x.txt", "y.txt", "fx.txt", "fy.txt"
        for name, values in zip(names, arrays):
            values = np.asarray(values)
            np.savetxt(work / name, values.reshape(len(values), -1), fmt="%.8g")
        command = [bcpd, "-xx.txt", "-yy.txt", "-Xfx.txt", "-Yfy.txt",
                   *options.split()]
        subprocess.run(command, cwd=work, check=True)
        result = np.loadtxt(work / "output_y.txt")
        if trajectory:
            subprocess.run([sys.executable, VIEWERS / "optpath.py",
                            work / ".optpath.bin"], check=True)
        return result
