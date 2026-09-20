#!/usr/bin/env python3
"""Run the three-level female or male GBCPD registration."""

import argparse

import numpy as np

from _common import compare_mesh, data_file, register

LEVELS = ((1.2, .5, 3), (.7, 0, .1), (.3, 0, .1))


if __name__ == "__main__":
    name = argparse.ArgumentParser(description=__doc__)
    name.add_argument("model", choices=("female", "male"))
    name = name.parse_args().model

    target = np.loadtxt(data_file(f"{name}-x.txt"))
    original = result = np.loadtxt(data_file(f"{name}-y.txt"))
    faces = data_file(f"{name}-triangles.txt")
    for beta, tau, gamma in LEVELS:
        result = register(target, result,
                          f"-Ggeodesic,{tau},{faces} -w0 -b{beta} -l100 "
                          f"-g{gamma} -J300 -K300 -p -ux "
                          "-c1e-6 -n500 -r1")
    compare_mesh(target, original, result, np.loadtxt(faces, dtype=int))
