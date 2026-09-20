#!/usr/bin/env python3
"""Run the grayscale or color image-registration DET demo."""

import argparse

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy.ndimage import median_filter

from _common import data_file, run


def read(name):
    image = np.asarray(Image.open(data_file(name)), dtype=float) / 255
    return image[..., None] if image.ndim == 2 else image


def fields(image):
    height, width = image.shape[:2]
    row, column = np.indices((height, width))
    points = np.column_stack((column.ravel() + 1, row.ravel() + 1))
    return points, image.reshape(-1, image.shape[2])


def warp(image, points, shape):
    height, width = shape[:2]
    output = np.zeros((height, width, image.shape[2]))
    xy = np.rint(points).astype(int) - 1
    valid = ((xy[:, 0] >= 0) & (xy[:, 0] < width) &
             (xy[:, 1] >= 0) & (xy[:, 1] < height))
    output[xy[valid, 1], xy[valid, 0]] = image.reshape(-1, image.shape[2])[valid]
    return output


def show(images, titles, gray=False):
    fig, axes = plt.subplots(1, len(images), figsize=(4 * len(images), 4))
    for ax, image, title in zip(axes, images, titles):
        ax.imshow(np.squeeze(np.clip(image, 0, 1)), cmap="gray" if gray else None)
        ax.set_title(title)
        ax.axis("off")
    fig.tight_layout()
    plt.show()


def color():
    target, source = read("antonio1-s.jpg"), read("antonio2-s.jpg")
    x, fx = fields(target)
    y, fy = fields(source)
    result = run(x, y, fx, fy,
                 "-A -l50 -w.1 -b1 -g1 -n500 -j1 -DB,5000,.03,.05")
    images = (source, target, warp(source, result, target.shape))
    show(tuple(np.rot90(image, -1) for image in images),
         ("Source", "Target", "Registered"))


def gray():
    target, source = read("111.jpg"), read("222.jpg")
    x, fx = fields(target)
    y, fy = fields(source)
    level1 = run(x, y, fx, fy,
                 "-A -n500 -l20 -b.9 -w0 -c1e-6 -g3 -DB,10000,.01 -r1")
    level2 = run(x, level1, fx, fy,
                 "-A -n500 -l20 -b.8 -w0 -c1e-6 -g.1 -DB,30000,.01 -r1")
    image1 = median_filter(warp(source, level1, target.shape)[..., 0], 5)
    image2 = median_filter(warp(source, level2, target.shape)[..., 0], 5)
    target, source = target[..., 0], source[..., 0]
    show((source, target, image1, image2, abs(source - target),
          abs(image1 - target), abs(image2 - target)),
         ("Source", "Target", "Level 1", "Level 2", "Before error",
          "Level 1 error", "Level 2 error"), True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("demo", choices=("gray", "color"))
    globals()[parser.parse_args().demo]()
