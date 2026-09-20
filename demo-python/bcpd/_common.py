"""Shared utilities for the BCPD Python demos."""

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
from _viewer import compare, show_cloud
from meshpath import compare as compare_mesh


def data_file(name):
    for folder in (ROOT / "data", ROOT / "data/local"):
        path = folder / name
        if path.exists():
            return path
    raise FileNotFoundError(f"Data file not found: {name}")


def execute(x_name, y_name, options, trajectory=False, faces=None):
    """Run BCPD in a temporary directory and visualize its result."""
    x, y = data_file(x_name), data_file(y_name)
    faces = data_file(faces) if faces else None
    options = options.format(faces=faces)
    bcpd = ROOT / ("win/bcpd.exe" if os.name == "nt" else "bcpd")
    if not bcpd.exists():
        raise FileNotFoundError("BCPD is not built; run make in the repository root")

    with tempfile.TemporaryDirectory(prefix="bcpd-demo-") as work:
        subprocess.run([bcpd, f"-x{x}", f"-y{y}", *options.split()],
                       cwd=work, check=True)
        if faces:
            files = ([Path(work) / ".optpath.bin"] if trajectory else
                     [x, y, Path(work) / "output_y.txt", faces])
            subprocess.run([sys.executable, VIEWERS / "meshpath.py", *files],
                           check=True)
        elif trajectory:
            subprocess.run([sys.executable, VIEWERS / "optpath.py",
                            Path(work) / ".optpath.bin"], check=True)
        else:
            compare(*(np.loadtxt(file) for file in
                      (x, y, Path(work) / "output_y.txt")))


def rigid_transform(x_name, y_name, options):
    """Run rigid BCPD and return its rotation and translation."""
    bcpd = ROOT / ("win/bcpd.exe" if os.name == "nt" else "bcpd")
    with tempfile.TemporaryDirectory(prefix="bcpd-demo-") as work:
        subprocess.run([bcpd, f"-x{data_file(x_name)}", f"-y{data_file(y_name)}",
                        *options.split()], cwd=work, check=True)
        return (np.loadtxt(Path(work) / "output_R.txt"),
                np.loadtxt(Path(work) / "output_t.txt"))


def register(target, source, options):
    """Register two arrays and return the transformed source."""
    bcpd = ROOT / ("win/bcpd.exe" if os.name == "nt" else "bcpd")
    with tempfile.TemporaryDirectory(prefix="bcpd-demo-") as work:
        work = Path(work)
        np.savetxt(work / "x.txt", target)
        np.savetxt(work / "y.txt", source)
        subprocess.run([bcpd, "-xx.txt", "-yy.txt", *options.split()],
                       cwd=work, check=True)
        return np.loadtxt(work / "output_y.txt")
