#!/usr/bin/env python3
"""Plot one or more GRABNES Brillouin-zone debug files."""

import argparse
from pathlib import Path

from bz_io import closed_polygon, load_table, parse_bz_file


DEFAULT_INPUTS = [
    "brillouin_zones_debug_K1.dat",
    "brillouin_zones_debug_K2.dat",
    "brillouin_zones_debug_K3.dat",
]


def build_plot(inputs, g_vectors=None):
    import matplotlib.pyplot as plt

    fig, axis = plt.subplots(figsize=(10, 8))
    colors = ("tab:blue", "tab:red", "tab:green", "tab:purple")

    for index, input_path in enumerate(inputs):
        vertices, k_point, b1, b2 = parse_bz_file(input_path)
        color = colors[index % len(colors)]
        label = f"Layer {index + 1}"
        polygon = closed_polygon(vertices)
        axis.plot(polygon[:, 0], polygon[:, 1], color=color, label=f"{label} BZ")
        axis.fill(polygon[:, 0], polygon[:, 1], color=color, alpha=0.12)
        if k_point is not None:
            axis.scatter(*k_point, color=color, marker="*", s=90)
        if index == 0:
            axis.quiver(0, 0, b1[0], b1[1], angles="xy", scale_units="xy", scale=1)
            axis.quiver(0, 0, b2[0], b2[1], angles="xy", scale_units="xy", scale=1)

    vectors = load_table(g_vectors) if g_vectors else None
    if vectors is not None:
        axis.scatter(vectors[:, 0], vectors[:, 1], s=14, color="0.35", alpha=0.6,
                     label=f"G vectors ({len(vectors)})")

    axis.scatter(0, 0, color="black", s=35, label="Origin")
    axis.set(xlabel=r"$k_x$ ($\AA^{-1}$)", ylabel=r"$k_y$ ($\AA^{-1}$)",
             title="GRABNES Brillouin-zone debug output")
    axis.set_aspect("equal")
    axis.grid(alpha=0.25)
    axis.legend()
    fig.tight_layout()
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="*", type=Path, help="BZ debug data files")
    parser.add_argument("--g-vectors", type=Path, help="optional G-vector table")
    parser.add_argument("--output", type=Path, default=Path("brillouin_zones.png"))
    parser.add_argument("--show", action="store_true", help="also open an interactive window")
    args = parser.parse_args()

    if not args.show:
        import matplotlib

        matplotlib.use("Agg")
    inputs = args.inputs or [Path(value) for value in DEFAULT_INPUTS]
    figure = build_plot(inputs, args.g_vectors)
    figure.savefig(args.output, dpi=180)
    print(f"Wrote {args.output}")
    if args.show:
        import matplotlib.pyplot as plt

        plt.show()


if __name__ == "__main__":
    main()
