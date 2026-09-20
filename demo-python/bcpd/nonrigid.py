#!/usr/bin/env python3
"""Run one of the original BCPD nonrigid-registration demos."""

import argparse

from _common import execute

CASES = {
    "armadillo-a": ("armadillo-x.txt", "armadillo-y.txt", "-w0 -b2 -l20 -g10 -J300 -K70 -p -c1e-6 -n500 -r1 -sY"),
    "armadillo-b": ("armadillo-y.txt", "armadillo-x.txt", "-w0 -b2 -l20 -g10 -J300 -K70 -p -c1e-6 -n500 -r1 -sY"),
    "bunny-a": ("bunny-y.txt", "bunny-x-co.txt", "-w0.1 -b2 -l200 -g3 -J300 -K70 -p -c1e-6 -n500 -r1 -sA"),
    "bunny-b": ("bunny-x-co.txt", "bunny-y.txt", "-w0.1 -b2 -l100 -g10 -J300 -K70 -p -c1e-6 -n500 -r1 -sA"),
    "dragon-a": ("dragon-y.txt", "dragon-x.txt", "-w0 -b2.5 -l4 -g10 -J300 -K50 -p -e.25 -c1e-6 -n90 -r1 -sY"),
    "dragon-b": ("dragon-x.txt", "dragon-y.txt", "-w0 -b2.5 -l4 -g10 -J300 -K50 -p -e.25 -c1e-6 -n90 -r1 -sY"),
    "face-a": ("face-y.txt", "face-x.txt", "-w0 -b0.3 -l1e4 -g10 -J300 -K150 -p -c1e-6 -n500 -r2 -sA"),
    "face-b": ("face-x.txt", "face-y.txt", "-w0 -b0.3 -l1e4 -g10 -J300 -K150 -p -c1e-6 -n500 -r2 -sA"),
    "fish-a": ("fish-x.txt", "fish-y.txt", "-w0 -b2 -l2 -g3 -c1e-6 -n500 -sA"),
    "fish-b": ("fish-y.txt", "fish-x.txt", "-w0 -b2 -l2 -g3 -c1e-6 -n500 -sA"),
    "monkey-a": ("monkey-y.txt", "monkey-x.txt", "-w0.1 -b2 -l200 -g10 -J300 -K70 -p -c1e-6 -n500 -r1 -sA"),
    "monkey-b": ("monkey-x.txt", "monkey-y.txt", "-w0.1 -b2 -l200 -g10 -J300 -K70 -p -c1e-6 -n500 -r1 -sA"),
}

def run(name):
    execute(*CASES[name], trajectory=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("demo", choices=CASES)
    run(parser.parse_args().demo)
