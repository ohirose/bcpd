"""Shared utilities for the DLD Python demos."""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "demo-python/viewers"))
from _viewer import compare


def data_file(name):
    for folder in (ROOT / "data", ROOT / "data/local"):
        path = folder / name
        if path.exists():
            return path
    raise FileNotFoundError(name)


def model(shapes, rank):
    """Return the mean and DLD covariance representation."""
    shapes = np.asarray(shapes)
    mean = shapes.mean(0)
    _, values, vectors = np.linalg.svd(
        (shapes - mean).reshape(len(shapes), -1), full_matrices=False)
    return mean, np.vstack((values[:rank] ** 2 / len(shapes),
                            vectors[:rank].T))


def register(target, source, options, covariance=None):
    """Run BCPD/DLD on arrays and return the transformed source."""
    bcpd = ROOT / ("win/bcpd.exe" if os.name == "nt" else "bcpd")
    with tempfile.TemporaryDirectory(prefix="dld-demo-") as work:
        work = Path(work)
        np.savetxt(work / "x.txt", target)
        np.savetxt(work / "y.txt", source)
        command = [bcpd, "-xx.txt", "-yy.txt", *options.split()]
        if covariance is not None:
            np.savetxt(work / "c.txt", covariance)
            command.append("-Cc.txt")
        subprocess.run(command, cwd=work, check=True)
        return np.loadtxt(work / "output_y.txt")


def fit(target, training, rank, dld_lambda, refine, groups):
    mean, covariance = model(training, rank)
    result = register(target, mean,
                      f"-un -l{dld_lambda} -w1e-2 -g1 -n50", covariance)
    if refine:
        result = register(target, result, "-ux -l10 -b.3 -w0 -g1 -n50")
    compare(target, mean, result, groups)
