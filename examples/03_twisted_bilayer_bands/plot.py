#!/usr/bin/env python3
"""Plot a compact GRABNES twisted-bilayer band calculation."""

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
OUTPUT = HERE / "twisted_bilayer_bands.png"


def read_bands(path: Path) -> tuple[list[float], list[list[float]], float]:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    number_of_bands, number_of_spins, number_of_points = map(int, lines[3].split())
    values = [float(value) for line in lines[4:] for value in line.split()]
    row_width = 1 + number_of_bands * number_of_spins
    if len(values) != row_width * number_of_points:
        raise ValueError(f"{path}: inconsistent band-data size")

    x_values: list[float] = []
    energies: list[list[float]] = [[] for _ in range(row_width - 1)]
    for offset in range(0, len(values), row_width):
        row = values[offset : offset + row_width]
        x_values.append(row[0])
        for band, energy in zip(energies, row[1:]):
            band.append(energy)
    return x_values, energies, float(lines[1].split()[-1])


def main() -> None:
    source = CALCULATED if CALCULATED.exists() else REFERENCE
    x_values, energies, path_end = read_bands(source)
    ratios = [0.0, 1.702760 / 4.028774, 3.177394 / 4.028774, 1.0]
    symmetry_positions = [path_end * ratio for ratio in ratios]

    figure, axis = plt.subplots(figsize=(6.4, 4.6))
    for band in energies:
        axis.plot(x_values, band, color="#2457a6", linewidth=0.7)
    for position in symmetry_positions:
        axis.axvline(position, color="0.82", linewidth=0.8, zorder=0)
    axis.axhline(0.0, color="0.45", linewidth=0.8, linestyle="--")
    axis.set_xticks(symmetry_positions, [r"$K$", r"$\Gamma$", r"$M$", r"$K'$" ])
    axis.set_xlim(0.0, path_end)
    axis.set_ylim(-3.0, 3.0)
    axis.set_xlabel("Wave-vector path")
    axis.set_ylabel("Energy (eV)")
    axis.set_title("Commensurate twisted-bilayer graphene bands")
    figure.tight_layout()
    figure.savefig(OUTPUT, dpi=180)
    print(f"Read {source.name}; wrote {OUTPUT}")


if __name__ == "__main__":
    main()
