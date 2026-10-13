"""A synthetic Fish control with identities known BEFORE registration.

Requires NumPy. Run `python landmark_check.py --help` for the two steps.
This is an in-sample teaching check, not anatomical or held-out validation.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np


def coordinates(path):
    values = np.loadtxt(path, ndmin=2)
    if values.shape != (91, 2) or not np.isfinite(values).all():
        raise ValueError(f"Expected 91 finite 2D Fish points in {path}")
    return values


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(fish, folder):
    target = coordinates(fish)
    # A known invertible shear, followed by similarity motion.
    source = target.copy()
    source[:, 1] += 0.18 * np.sin(1.3 * target[:, 0])
    angle = np.deg2rad(10)
    rotation = np.array([[np.cos(angle), -np.sin(angle)],
                         [np.sin(angle), np.cos(angle)]])
    source = 1.04 * source @ rotation.T + np.array([0.22, -0.14])

    # Independent permutations: matching row numbers are NOT matching IDs.
    rng = np.random.default_rng(42)
    target_order = rng.permutation(len(target))
    source_order = rng.permutation(len(source))
    target_rows = np.argsort(target_order) + 1
    source_rows = np.argsort(source_order) + 1
    # Five fixed controls, chosen before fitting (zero-based original rows).
    controls = dict(zip([23, 58, 6, 36, 45], "ABCDE"))

    folder.mkdir(parents=True, exist_ok=False)
    np.savetxt(folder / "target.txt", target[target_order], fmt="%.17g")
    np.savetxt(folder / "source.txt", source[source_order], fmt="%.17g")
    with (folder / "identities.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["point_id", "label", "target_row", "source_row"])
        for i in range(len(target)):
            writer.writerow([f"P{i+1:03d}", controls.get(i, ""),
                             target_rows[i], source_rows[i]])
    metadata = {
        "kind": "synthetic, in-sample; identities withheld from registration",
        "units": "original Fish coordinate units, not verified millimetres",
        "warp": "1.04 * (x, y + 0.18*sin(1.3*x)) @ R(10deg).T + (0.22,-0.14)",
        "permutation_seed": 42,
        "fish_sha256": digest(fish),
        "sha256": {name: digest(folder / name) for name in
                   ["target.txt", "source.txt", "identities.csv"]},
    }
    (folder / "control.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Prepared {folder}: 91 points and an evaluation-only identity map.")
    print("Run BCPD on target.txt and source.txt, NOT on identities.csv.")


def measure(folder, result):
    metadata = json.loads((folder / "control.json").read_text())
    for name, expected in metadata["sha256"].items():
        if digest(folder / name) != expected:
            raise ValueError(f"Control file changed: {name}; prepare a fresh control")
    target = coordinates(folder / "target.txt")
    source = coordinates(folder / "source.txt")
    transformed = coordinates(result)
    with (folder / "identities.csv").open(newline="") as stream:
        identities = list(csv.DictReader(stream))
    if len({row["point_id"] for row in identities}) != 91:
        raise ValueError("Expected 91 unique control IDs")
    ti = np.array([int(row["target_row"]) - 1 for row in identities])
    si = np.array([int(row["source_row"]) - 1 for row in identities])
    if not (np.array_equal(np.sort(ti), np.arange(91)) and
            np.array_equal(np.sort(si), np.arange(91))):
        raise ValueError("Identity map must contain each source and target row once")
    before = np.linalg.norm(source[si] - target[ti], axis=1)
    after = np.linalg.norm(transformed[si] - target[ti], axis=1)
    print("Synthetic control; original coordinate units; in-sample, NOT held-out.")
    print("Label ID    Target row Source row       Before        After")
    selected = [(i, row) for i, row in enumerate(identities) if row["label"]]
    selected.sort(key=lambda item: item[1]["label"])
    if [row["label"] for _, row in selected] != list("ABCDE"):
        raise ValueError("Expected exactly the five labels A-E")
    for i, row in selected:
        print(f"{row['label']:5} {row['point_id']:5} {ti[i]+1:10d} {si[i]+1:10d} "
              f"{before[i]:12.6g} {after[i]:12.6g}")
    indices = [i for i, _ in selected]
    for label, values in [("Five-control mean", indices), ("All 91 mean", slice(None))]:
        print(f"{label:38} {np.mean(before[values]):12.6g} {np.mean(after[values]):12.6g}")
    print(f"{'All 91 RMS':38} {np.sqrt(np.mean(before**2)):12.6g} "
          f"{np.sqrt(np.mean(after**2)):12.6g}")
    print(f"{'All 91 maximum':38} {np.max(before):12.6g} {np.max(after):12.6g}")
    print("Errors near zero are limited by the precision of the saved coordinates.")
    print("Use this control's final y file without reordering its rows.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    create = commands.add_parser("prepare", help="make a NEW synthetic control folder")
    create.add_argument("fish", type=Path, help="distributed fish-x.txt")
    create.add_argument("folder", type=Path, help="new output directory; must not exist")
    check = commands.add_parser("measure", help="measure final y using the saved identity map")
    check.add_argument("folder", type=Path)
    check.add_argument("result", type=Path, help="final denormalized BCPD *_y.txt")
    args = parser.parse_args()
    try:
        if args.command == "prepare":
            prepare(args.fish, args.folder)
        else:
            measure(args.folder, args.result)
    except (OSError, ValueError, KeyError) as error:
        parser.exit(1, f"Error: {error}\n")


if __name__ == "__main__":
    main()
