#!/usr/bin/env python3
"""Plot the calculated or reference twisted-bilayer density of states."""

from pathlib import Path

try:
    import matplotlib.pyplot as plt
except ImportError as exc:
    raise SystemExit(
        "plot.py requires matplotlib. Install it with: python3 -m pip install matplotlib"
    ) from exc


HERE = Path(__file__).resolve().parent
CALCULATED = HERE / "generate.diag.DOS"
REFERENCE = HERE / "reference" / "dos.dat"
OUTPUT = HERE / "twisted_bilayer_dos.png"


def read_dos(path: Path) -> tuple[list[float], list[float]]:
    energy: list[float] = []
    density: list[float] = []
    for line_number, line in enumerate(path.read_text().splitlines(), start=1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        values = line.split()
        if len(values) < 2:
            raise ValueError(f"{path}:{line_number}: expected at least two columns")
        energy.append(float(values[0]))
        density.append(sum(float(value) for value in values[1:]))
    return energy, density


def main() -> None:
    source = CALCULATED if CALCULATED.exists() else REFERENCE
    energy, density = read_dos(source)
    figure, axis = plt.subplots(figsize=(6.4, 4.2))
    axis.plot(energy, density, color="#7b3294", linewidth=1.6)
    axis.axvline(0.0, color="0.5", linestyle="--", linewidth=0.8)
    axis.set_xlabel("Energy (eV)")
    axis.set_ylabel("Density of states (states/eV)")
    axis.set_title("Commensurate twisted-bilayer graphene DOS")
    figure.tight_layout()
    figure.savefig(OUTPUT, dpi=180)
    print(f"Read {source.name}; wrote {OUTPUT}")


if __name__ == "__main__":
    main()
