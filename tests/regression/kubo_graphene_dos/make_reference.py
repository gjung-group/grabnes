#!/usr/bin/env python3
"""Exact density of states of nearest-neighbor graphene for the Kubo check.

The reference is analytic, not solver output: the two bands are +-t|f(k)| with
f(k) = 1 + exp(i k.a1) + exp(i k.a2), and the recursion (continued-fraction)
DOS of GRABNES is a Lorentzian-broadened DOS with half width Epsilon. The
solver works in units of the hopping t (= TB.Hopping), so energies, Epsilon and
the DOS (states per atom and per t) are expressed in that unit.

Requires NumPy. Run it only to regenerate reference/exact_dos.dat after
changing Gendata.in.
"""

from pathlib import Path

import numpy as np

HOPPING = 3.1          # TB.Hopping (eV)
EPSILON = 0.05         # Epsilon (eV)
E_MIN, E_MAX = -9.0, 9.0
N_ENERGY = 201         # NumberofEnergyPoints; the solver writes N_ENERGY + 1 rows
N_K = 1500             # k-points per direction for the Brillouin-zone average

eta = EPSILON / HOPPING
energies = (E_MIN + (E_MAX - E_MIN) * np.arange(N_ENERGY + 1) / N_ENERGY) / HOPPING

u = (np.arange(N_K) + 0.5) / N_K
k1, k2 = np.meshgrid(u, u)
f = np.abs(1.0 + np.exp(2j * np.pi * k1) + np.exp(2j * np.pi * k2)).ravel()


def lorentzian(x):
    return eta / np.pi / (x * x + eta * eta)


dos = np.array([0.5 * (lorentzian(e - f).mean() + lorentzian(e + f).mean()) for e in energies])

out = Path(__file__).resolve().parent / "reference" / "exact_dos.dat"
with open(out, "w") as handle:
    for e, d in zip(energies, dos):
        handle.write(f"{e:24.16e} {d:24.16e}\n")
print(f"wrote {out}")
