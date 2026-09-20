#!/usr/bin/env python3
"""Display a BCPD point-set trajectory."""

import sys

import numpy as np

from _viewer import compare, trajectory


def load(filename):
    with open(filename, "rb") as file:
        n, d, m, steps = map(int, np.fromfile(file, "<i4", 4))
        path = np.fromfile(file, "<f8", d * m * steps)
        target = np.fromfile(file, "<f8", d * n)
    return (target.reshape(d, n, order="F").T,
            path.reshape(d, m, steps, order="F").transpose(2, 1, 0))


if __name__ == "__main__":
    target, path = load(sys.argv[1] if len(sys.argv) > 1 else ".optpath.bin")
    trajectory(target, path)
    compare(target, path[0], path[-1])
