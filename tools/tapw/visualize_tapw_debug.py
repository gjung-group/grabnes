#!/usr/bin/env python3
"""Create a compact overview of TAPW reciprocal-space debug output."""

import argparse
from pathlib import Path

from bz_io import closed_polygon, load_table, parse_bz_file


def plot_debug_directory(directory):
    import matplotlib.pyplot as plt

    directory = Path(directory)
    bz_files = sorted(directory.glob("brillouin_zones_debug_K*.dat"))
    if not bz_files:
        raise FileNotFoundError(f"no brillouin_zones_debug_K*.dat files in {directory}")

    g_vectors = load_table(directory / "g_vectors_debug.dat")
    k_path = load_table(directory / "kpath_debug.dat_absolute")
    fig, axes = plt.subplots(1, 2, figsize=(13, 6))
    reciprocal_axis, path_axis = axes

    for index, path in enumerate(bz_files):
        vertices, k_point, _, _ = parse_bz_file(path)
        polygon = closed_polygon(vertices)
        reciprocal_axis.plot(polygon[:, 0], polygon[:, 1], label=f"Layer {index + 1}")
        if k_point is not None:
            reciprocal_axis.scatter(*k_point, marker="*", s=80)

    if g_vectors is not None:
        reciprocal_axis.scatter(g_vectors[:, 0], g_vectors[:, 1], s=12, color="0.4",
                                alpha=0.55, label=f"G vectors ({len(g_vectors)})")
    reciprocal_axis.set_title("Brillouin zones and TAPW basis")
    reciprocal_axis.set_aspect("equal")
    reciprocal_axis.grid(alpha=0.25)
    reciprocal_axis.legend()

    if k_path is None:
        path_axis.text(0.5, 0.5, "No kpath_debug.dat_absolute", ha="center", va="center")
    else:
        path_axis.plot(k_path[:, 0], k_path[:, 1], marker="o", markersize=3)
        path_axis.set_aspect("equal")
    path_axis.set_title("Absolute k path")
    path_axis.grid(alpha=0.25)

    for axis in axes:
        axis.set_xlabel(r"$k_x$ ($\AA^{-1}$)")
        axis.set_ylabel(r"$k_y$ ($\AA^{-1}$)")
    fig.tight_layout()
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", nargs="?", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path, default=Path("tapw_debug.png"))
    parser.add_argument("--show", action="store_true", help="also open an interactive window")
    args = parser.parse_args()

    if not args.show:
        import matplotlib

        matplotlib.use("Agg")
    figure = plot_debug_directory(args.directory)
    figure.savefig(args.output, dpi=180)
    print(f"Wrote {args.output}")
    if args.show:
        import matplotlib.pyplot as plt

        plt.show()


if __name__ == "__main__":
    main()
