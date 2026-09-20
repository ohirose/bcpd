#!/usr/bin/env python3
"""Run one of the original rigid BCPD scan-pair demos."""

import argparse

from _common import execute

SERIES = {
    "apartment": ("apartment_", 0, "0.1", "B,2000,0.08"),
    "chef": ("chef_view", 0, "0.2", "B,10000,0.04"),
    "para": ("parasaurolophus_view", 0, "0.2", "B,10000,0.04"),
    "stairs": ("stairs_", 1, "0.1", "B,2000,0.08"),
    "trex": ("T-rex_view", 0, "0.2", "B,10000,0.04"),
}
CASES = [f"{name}-{chr(97 + i)}" for name in SERIES for i in range(9)]


def preset(name):
    series, letter = name.rsplit("-", 1)
    prefix, start, omega, downsample = SERIES[series]
    i = ord(letter) - ord("a")
    x_index = start + i + 1
    y_index = start + i if series in ("apartment", "stairs") else x_index + 1
    gamma, extra = "1", ""

    if name == "apartment-d":
        gamma, extra = "0.1", "-e0.3 -un "
    elif name == "apartment-i":
        omega, gamma, downsample = "0.05", "3", "B,5000,0.08"
    elif name == "trex-a":
        downsample = "B,10000,0.08"
    elif name == "trex-d":
        omega, downsample = "0.3", "B,10000,0.08"
    elif name == "trex-g":
        omega = "0.3"

    x = f"{prefix}{x_index:03d}.txt"
    y = f"{prefix}{y_index:03d}.txt"
    downsample = "" if name == "para-c" else f"-D{downsample} "
    options = (f"-w{omega} -g{gamma} -J300 -p {extra}-f0.3 {downsample}"
               "-c1e-6 -n500 -r1 -Tr -sT")
    return x, y, options


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("demo", choices=CASES)
    execute(*preset(parser.parse_args().demo))
