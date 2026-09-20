#!/usr/bin/env python3
"""Fit a leave-one-out DLD hand model."""

import argparse

import numpy as np

from _common import data_file, fit


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("hand", type=int, choices=range(1, 41), nargs="?", default=1)
    parser.add_argument("--refine", action="store_true", help="fine-tune with BCPD")
    args = parser.parse_args()
    load = lambda hand: np.loadtxt(data_file(f"hand{hand:03d}.txt"))
    training = [load(hand) for hand in range(1, 41) if hand != args.hand]
    fit(load(args.hand), training, 38, .1, args.refine, (range(56),))
