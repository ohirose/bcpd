#!/usr/bin/env python3
"""Fit a leave-one-person-out DLD face model."""

import argparse

import numpy as np

from _common import data_file, fit

GROUPS = (range(13), (*range(13, 21), 13), (*range(21, 29), 21),
          range(29, 34), range(34, 39), (*range(39, 47), 39), range(47, 58))


def face(person, expression):
    for sex in "mf":
        try:
            return np.loadtxt(data_file(
                f"{person:02d}-{expression}{sex}-opa.txt"))
        except FileNotFoundError:
            pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("person", type=int, choices=range(1, 41), nargs="?",
                        default=1)
    parser.add_argument("expression", type=int, choices=range(1, 7), nargs="?",
                        default=1)
    parser.add_argument("--refine", action="store_true", help="fine-tune with BCPD")
    args = parser.parse_args()
    training = [face(person, expression) for person in range(1, 41)
                if person != args.person for expression in range(1, 7)]
    fit(face(args.person, args.expression), training, 50, .5 if args.refine else 1,
        args.refine, GROUPS)
