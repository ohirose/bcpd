#!/usr/bin/env python3
"""Run one of the original BCPD++ demos."""

import argparse

from _common import execute

BASE = "-w0 -b2 -l50 -g10 -J300 -K70 -p -ue -c1e-6 -n500 -r1"
CASES = {
    "armadillo": ("armadillo-x.txt", "armadillo-y.txt", f"{BASE} -f0.3 -DB,5000,0.02"),
    "asian-dragon": ("asiandragon-x.txt", "asiandragon-y.txt", f"{BASE} -DB,50000,0.08"),
    "dragon": ("dragon-x.txt", "dragon-y.txt", f"{BASE} -DB,10000,0.08"),
    "lucy": ("lucy-x.txt", "lucy-y.txt", f"{BASE} -DB,50000,0.08"),
}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("demo", choices=CASES)
    execute(*CASES[parser.parse_args().demo])
