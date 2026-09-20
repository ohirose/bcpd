#!/usr/bin/env python3
"""Run a DET shape-localization or shape-registration demo."""

import argparse

import numpy as np

from _common import compare, data_file, run

LOCALIZATION = {
    "localization-1": ("y.txt", "skull.txt", "fy2.txt", "fy2.skull.txt",
                       "-ux -Tr -w.9 -j.3 -g1 -n500 -sY -c1e-7 -r1"),
    "localization-2": ("skull.txt", "y.txt", "fy2.skull.txt", "fy2.txt",
                       "-ux -Tr -w.1 -g3 -n500 -sY -r1"),
}


def load(name):
    return np.loadtxt(data_file(name))


def localization(name):
    x, y, fx, fy, options = LOCALIZATION[name]
    run(load(x), load(y), load(fx), load(fy), options, trajectory=True)


def registration():
    target, source = load("tr_reg_003.vert.txt"), load("tr_reg_001.vert.txt")
    result = run(target, source, load("tr_reg_003.wks.c30.txt"),
                 load("tr_reg_001.wks.c30.txt"),
                 f"-Ggeo,.1,{data_file('tr_reg.face.txt')} -ux -n200 -c1e-6 "
                 "-A -e.4 -l100 -b1 -w.1 -g1 -r1 -sY -DB,2000,.05")
    compare(target, source, result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("demo", choices=(*LOCALIZATION, "registration"))
    demo = parser.parse_args().demo
    localization(demo) if demo in LOCALIZATION else registration()
