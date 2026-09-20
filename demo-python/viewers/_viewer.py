"""Matplotlib for 2D points; Open3D for 3D points."""

import time

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

BLUE, RED, BLACK = (0, .45, .74), (.85, .1, .1), (0, 0, 0)


def open3d():
    try:
        import open3d as o3d
        return o3d
    except ImportError:
        raise SystemExit("Open3D is not installed; run: python -m pip install open3d")


def cloud(vertices, color):
    o3d = open3d()
    result = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(vertices))
    result.paint_uniform_color(color)
    return result


def visualizer(title, *geometry, left=50, width=800):
    o3d = open3d()
    viewer = o3d.visualization.Visualizer()
    viewer.create_window(title, width=width, height=700, left=left, top=50)
    for item in geometry:
        viewer.add_geometry(item)
    viewer.get_render_option().background_color = np.ones(3)
    viewer.get_render_option().point_size = 2
    return viewer


def show_views(*views):
    """Show independent views; closing any window closes them all."""
    width = 600 if len(views) > 2 else 800
    viewers = [visualizer(title, *geometry, left=20 + i * (width + 20),
                          width=width)
               for i, (title, geometry) in enumerate(views)]
    running = True
    while running:
        for viewer in viewers:
            running = viewer.poll_events()
            if not running:
                break
            viewer.update_renderer()
        time.sleep(.01)
    for viewer in viewers:
        viewer.destroy_window()


def show_pair(before, after):
    show_views(("Before Registration", before), ("After Registration", after))


def show_cloud(points, colors):
    """Display one colored point cloud."""
    o3d = open3d()
    geometry = cloud(points, BLACK)
    geometry.colors = o3d.utility.Vector3dVector(colors)
    o3d.visualization.draw_geometries([geometry], window_name="3D Reconstruction")


def compare(target, source, result, groups=None):
    """Show registration before and after, choosing the viewer by dimension."""
    if target.shape[1] == 3:
        show_pair((cloud(target, BLUE), cloud(source, RED)),
                  (cloud(target, BLUE), cloud(result, RED)))
        return
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    for ax, points, title in zip(axes, (source, result),
                                 ("Before Registration", "After Registration")):
        if groups:
            for group in groups:
                ax.plot(*target[list(group)].T, "-o", ms=3, color="tab:blue")
                ax.plot(*points[list(group)].T, "-o", ms=3, color="tab:red")
        else:
            ax.plot(*target.T, ".", ms=2, label="Target")
            ax.plot(*points.T, ".", ms=2, color="tab:red", label="Source")
        ax.set(title=title, aspect="equal")
    if not groups:
        axes[0].legend()
    fig.tight_layout()
    plt.show()


def animate3d(target, path, moving):
    viewer = visualizer("Optimization Trajectory", cloud(target, BLUE), moving)
    frame, previous = 1, time.monotonic()

    def update(viewer):
        nonlocal frame, previous
        now = time.monotonic()
        if frame < len(path) and now - previous >= .03:
            moving.points = open3d().utility.Vector3dVector(path[frame])
            viewer.update_geometry(moving)
            frame, previous = frame + 1, now
        return False

    viewer.register_animation_callback(update)
    viewer.run()
    viewer.destroy_window()


def trajectory(target, path):
    """Animate a 2D path with Matplotlib or a 3D path with Open3D."""
    if target.shape[1] == 3:
        animate3d(target, path, cloud(path[0], RED))
        return

    low = np.minimum(target.min(0), path.min((0, 1)))
    high = np.maximum(target.max(0), path.max((0, 1)))
    pad = np.maximum(.05 * (high - low), 1e-9)
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(*target.T, ".", ms=2, label="Target")
    moving, = ax.plot(*path[0].T, ".", ms=2, color="tab:red", label="Source")
    ax.set(xlim=(low[0] - pad[0], high[0] + pad[0]),
           ylim=(low[1] - pad[1], high[1] + pad[1]), aspect="equal")
    ax.legend()

    def update(frame):
        moving.set_data(*path[frame].T)
        ax.set_title(f"Optimization trajectory ({frame + 1}/{len(path)})")
        return moving,

    animation = FuncAnimation(fig, update, len(path), interval=10, repeat=False)
    plt.show()
