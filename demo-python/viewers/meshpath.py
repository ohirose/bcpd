#!/usr/bin/env python3
"""Display a GBCPD mesh trajectory or comparison with Open3D."""

import sys

import numpy as np

from _viewer import BLACK, BLUE, animate3d, cloud, open3d, show_pair, show_views


def load(filename):
    with open(filename, "rb") as file:
        n, d, m, steps = map(int, np.fromfile(file, "<i4", 4))
        path = np.fromfile(file, "<f8", d * m * steps)
        target = np.fromfile(file, "<f8", d * n)
        rows, columns = map(int, np.fromfile(file, "<i4", 2))
        faces = np.fromfile(file, "<i4", rows * columns)
    return (target.reshape(d, n, order="F").T,
            path.reshape(d, m, steps, order="F").transpose(2, 1, 0),
            faces.reshape(columns, rows, order="F").T)


def wire(vertices, faces):
    o3d = open3d()
    faces = np.asarray(faces, dtype=int).reshape(-1, 3).copy()
    if faces.min() == 1:
        faces -= 1
    mesh = o3d.geometry.TriangleMesh(o3d.utility.Vector3dVector(vertices),
                                     o3d.utility.Vector3iVector(faces))
    result = o3d.geometry.LineSet.create_from_triangle_mesh(mesh)
    result.paint_uniform_color(BLACK)
    return result


def surface(vertices, faces):
    o3d = open3d()
    mesh = o3d.geometry.TriangleMesh(o3d.utility.Vector3dVector(vertices),
                                     o3d.utility.Vector3iVector(faces))
    mesh.compute_vertex_normals()
    mesh.paint_uniform_color((.7, .7, .7))
    return mesh


def transfer(target, target_faces, source, source_faces, result):
    show_views(("Source", (surface(source, source_faces),)),
               ("Target", (surface(target, target_faces),)),
               ("Deformed Source", (surface(result, source_faces),)))


def compare(target, before, after, faces):
    show_pair((cloud(target, BLUE), wire(before, faces)),
              (cloud(target, BLUE), wire(after, faces)))


def main(files):
    if len(files) == 1:
        target, path, faces = load(files[0])
        animate3d(target, path, wire(path[0], faces))
        compare(target, path[0], path[-1], faces)
    elif len(files) == 4:
        target, before, after = (np.loadtxt(file) for file in files[:3])
        faces = np.loadtxt(files[3], dtype=int)
        compare(target, before, after, faces)
    else:
        raise SystemExit("usage: meshpath.py OPTPATH | X Y REGISTERED FACES")


if __name__ == "__main__":
    main(sys.argv[1:])
