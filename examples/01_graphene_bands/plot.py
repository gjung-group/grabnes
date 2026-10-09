#!/usr/bin/env python3
"""Plot a GRABNES band file from this example."""

from pathlib import Path

try:
    import matplotlib.pyplot as plt
except ImportError as exc:
    raise SystemExit(
        "plot.py requires matplotlib. Install it with: python3 -m pip install matplotlib"
    ) from exc


HERE = Path(__file__).resolve().parent
CALCULATED = HERE / "generate.bands"
REFERENCE = HERE / "reference" / "bands.dat"
OUTPUT = HERE / "graphene_bands.png"


def read_bands(path: Path) -> tuple[list[float], list[list[float]]]:
    """Read the compact GRABNES .bands format."""
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if len(lines) < 5:
        raise ValueError(f"{path} does not contain a complete GRABNES band data set")

    metadata = lines[3].split()
    if not metadata:
        raise ValueError(f"Could not read band metadata from {path}")
    number_of_bands = int(metadata[0])

    x_values: list[float] = []
    energies: list[list[float]] = [[] for _ in range(number_of_bands)]
    for line_number, line in enumerate(lines[4:], start=5):
        values = [float(value) for value in line.split()]
        if len(values) != number_of_bands + 1:
            raise ValueError(
                f"{path}:{line_number}: expected {number_of_bands + 1} columns, "
                f"found {len(values)}"
            )
        x_values.append(values[0])
        for band, energy in zip(energies, values[1:]):
            band.append(energy)

    return x_values, energies


def main() -> None:
    source = CALCULATED if CALCULATED.exists() else REFERENCE
    x_values, energies = read_bands(source)

    figure, axis = plt.subplots(figsize=(6.4, 4.2))
    for band in energies:
        axis.plot(x_values, band, color="#2457a6", linewidth=1.8)

    symmetry_positions = [0.0, 1.702760, 3.177394, 4.028774]
    for position in symmetry_positions:
        axis.axvline(position, color="0.82", linewidth=0.8, zorder=0)
    axis.axhline(0.0, color="0.45", linewidth=0.8, linestyle="--")
    axis.set_xticks(symmetry_positions, [r"$K$", r"$\Gamma$", r"$M$", r"$K'$" ])
    axis.set_xlim(symmetry_positions[0], symmetry_positions[-1])
    axis.set_xlabel("Wave-vector path")
    axis.set_ylabel("Energy (eV)")
    axis.set_title("Pristine graphene band structure")
    figure.tight_layout()
    figure.savefig(OUTPUT, dpi=180)
    print(f"Read {source.name}; wrote {OUTPUT}")


if __name__ == "__main__":
    main()
