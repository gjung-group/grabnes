#!/usr/bin/env python3
"""Hermann four-index commensurate cells of two hexagonal lattices, in exact integer arithmetic.

Convention (the one of GRABNES, ``Structure.MoireCellParameters m1 m2 m3 m4``):

* the lattice vectors of a layer are a1 = a (1, 0) and a2 = a (1/2, sqrt(3)/2) of that layer, rotated with the layer;
* the supercell vector is A1 = m1 a1 + m2 a2 on the first lattice and the same vector is m3 a1' + m4 a2' on the
  second; A2 is A1 turned by +60 degrees.

With tau = exp(i pi/3) a lattice vector m a1 + n a2 is the Eisenstein integer m + n tau (times a). The pair is
alpha = m1 + m2 tau, beta = m3 + m4 tau. The second lattice is then (a alpha/beta) Z[tau]: it is turned by
theta = arg(alpha/beta) (counterclockwise, signed) and scaled by |alpha|/|beta| with respect to the first.

Facts used (Z[tau] = Z[omega] is a Euclidean domain, units = the six powers of tau):

* every pair of nonzero alpha, beta is a valid commensurate cell; alpha = 0 or beta = 0 is degenerate;
* the coincidence lattice of the two layers is alpha' Z[tau] with alpha' = alpha/gcd(alpha, beta): the cell is
  PRIMITIVE iff gcd(alpha, beta) is a unit, and it holds N(gcd) primitive cells (N = norm = |.|^2);
* (u alpha, u beta) for a unit u is the same cell (A1 turned by a multiple of 60 degrees); exchanging the layers
  gives (beta, alpha) and -theta; complex conjugation (a mirror) gives (conj alpha, conj beta) and -theta.

Run ``python3 hermann.py --selftest`` for the checks listed in the docstring of ``selftest``.
"""

import argparse
import math
from fractions import Fraction

TAU = complex(0.5, math.sqrt(3.0) / 2.0)


# ---------------------------------------------------------------- Eisenstein integers as (m, n) = m + n tau
def mul(x, y):
    """(a + b tau)(c + d tau), with tau^2 = tau - 1."""
    a, b = x
    c, d = y
    return (a * c - b * d, a * d + b * c + b * d)


def conj(x):
    """conj(a + b tau) = a + b (1 - tau) = (a + b) - b tau."""
    a, b = x
    return (a + b, -b)


def norm(x):
    a, b = x
    return a * a + a * b + b * b


def divmod_e(x, y):
    """Euclidean division x = q y + r with N(r) < N(y)."""
    num = mul(x, conj(y))
    n = norm(y)
    q = (round(Fraction(num[0], n)), round(Fraction(num[1], n)))
    best = None
    for da in (-1, 0, 1):          # nearest lattice point in the skew basis: check the neighbours
        for db in (-1, 0, 1):
            qq = (q[0] + da, q[1] + db)
            qy = mul(qq, y)
            r = (x[0] - qy[0], x[1] - qy[1])
            if best is None or norm(r) < norm(best[1]):
                best = (qq, r)
    assert norm(best[1]) < n
    return best


def exact_div(x, y):
    q, r = divmod_e(x, y)
    assert r == (0, 0), (x, y)
    return q


def gcd_e(x, y):
    while y != (0, 0):
        _, r = divmod_e(x, y)
        x, y = y, r
    return x


UNITS = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]   # tau^k, k = 0..5


def canonical(x):
    """The unit multiple of x with argument in [0, 60) degrees: m > 0, n >= 0."""
    for u in UNITS:
        y = mul(u, x)
        if y[0] > 0 and y[1] >= 0:
            return y, u
    raise ValueError(x)


# ---------------------------------------------------------------- cells
def to_complex(x):
    return x[0] + x[1] * TAU


def describe(m):
    """Everything about the cell (m1, m2, m3, m4): angle, ratio, primitivity, primitive indices."""
    al, be = (m[0], m[1]), (m[2], m[3])
    if norm(al) == 0 or norm(be) == 0:
        raise ValueError('degenerate: (m1, m2) and (m3, m4) must both be nonzero')
    g = gcd_e(al, be)
    alp, bep = exact_div(al, g), exact_div(be, g)
    lam = to_complex(al) / to_complex(be)
    theta = math.degrees(math.atan2(lam.imag, lam.real))
    return dict(theta=theta, ratio=abs(lam), n1=norm(al), n2=norm(be), mult=norm(g),
                primitive=(norm(g) == 1), prim=(alp[0], alp[1], bep[0], bep[1]),
                prim_n1=norm(alp), prim_n2=norm(bep))


def same_cell(m, mp):
    """True when two index sets describe the same structure: (alpha, beta) = u (alpha', beta')."""
    for u in UNITS:
        if mul(u, (mp[0], mp[1])) == (m[0], m[1]) and mul(u, (mp[2], mp[3])) == (m[2], m[3]):
            return True
    return False


def legacy_angle(m):
    """The unsigned angle printed by CellGet (acos formula, kept for backward compatibility)."""
    g = m[0] ** 2 + m[1] ** 2 + m[0] * m[1]
    d = math.sqrt((m[2] ** 2 + m[3] ** 2 + m[2] * m[3]) / g)
    c = (2 * m[0] * m[2] + 2 * m[1] * m[3] + m[0] * m[3] + m[1] * m[2]) / (2 * d * g)
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


# ---------------------------------------------------------------- searches
def search_tbg_hermann(nmax_cells):
    """All primitive commensurate cells of two IDENTICAL lattices with N(alpha') <= nmax_cells, from the four indices.

    Identical lattices: |alpha| = |beta|, i.e. N(alpha) = N(beta). Enumerate beta in the 60-degree sector (one
    representative per unit class) and every alpha of the same norm, keep the coprime pairs. Complete: every
    primitive pair is (u alpha, u beta) of exactly one pair met here.
    """
    by_norm = {}
    r = int(math.isqrt(4 * nmax_cells // 3)) + 2
    for a in range(-r, r + 1):
        for b in range(-r, r + 1):
            n = norm((a, b))
            if 0 < n <= nmax_cells:
                by_norm.setdefault(n, []).append((a, b))
    cells = []
    for n, vecs in by_norm.items():
        betas = [v for v in vecs if v[0] > 0 and v[1] >= 0]
        for be in betas:
            for al in vecs:
                if norm(gcd_e(al, be)) == 1:
                    cells.append((al[0], al[1], be[0], be[1]))
    return cells


def search_hermann(lam_target, tol, nmax_cells):
    """Commensurate primitive cells with |alpha/beta - lam_target| <= tol |lam_target| and N(alpha) <= nmax_cells.

    General lattice ratio and angle (heterostructures). Enumerate beta in the 60-degree sector with
    N(beta) <= nmax_cells / |lam|^2 (+1), and all alpha in the disc |alpha - lam beta| <= tol |lam| |beta|.
    Returned sorted by cell size (smallest first), then by the error.
    """
    lam = complex(lam_target)
    out = []
    bmax = int(math.isqrt(int(4 * nmax_cells / abs(lam) ** 2 / 3) + 4)) + 2
    for m3 in range(1, bmax + 1):
        for m4 in range(0, bmax + 1):
            be = (m3, m4)
            if norm(be) * abs(lam) ** 2 > nmax_cells * (1 + tol) ** 2 + 1:
                continue
            z = lam * to_complex(be)
            rad = tol * abs(z)
            # alpha = a + b tau near z: b = Im z / (sqrt3/2), a = Re z - b/2
            b0 = z.imag / (math.sqrt(3) / 2)
            a0 = z.real - b0 / 2
            span = int(rad * 2 / math.sqrt(3)) + 2
            for a in range(int(math.floor(a0)) - span, int(math.ceil(a0)) + span + 1):
                for b in range(int(math.floor(b0)) - span, int(math.ceil(b0)) + span + 1):
                    al = (a, b)
                    if norm(al) == 0 or norm(al) > nmax_cells:
                        continue
                    if abs(to_complex(al) - z) > rad:
                        continue
                    if norm(gcd_e(al, be)) != 1:
                        continue
                    lamc = to_complex(al) / to_complex(be)
                    out.append((norm(al), abs(lamc - lam) / abs(lam), (a, b, m3, m4)))
    out.sort()
    return out


def mele_angles(nmax_atoms):
    """Mele, PRB 81, 161405 (2010): cos(theta) = (n^2 + 4 n m + m^2) / (2 (n^2 + n m + m^2)), coprime m > n >= 1
    (plus theta = 0 is excluded). Atoms of the primitive cell: 4 (n^2 + n m + m^2) when (m - n) is not divisible by
    3, 4 (n^2 + n m + m^2)/3 when it is (the SE-even / SE-odd distinction)."""
    out = {}
    for m in range(1, nmax_atoms):
        if 4 * (m * m) / 3 > nmax_atoms * 3:
            break
        for n in range(1, m):
            if math.gcd(m, n) != 1:
                continue
            s = n * n + n * m + m * m
            nat = 4 * s // 3 if (m - n) % 3 == 0 else 4 * s
            if nat > nmax_atoms:
                continue
            c = (n * n + 4 * n * m + m * m) / (2 * s)
            th = math.degrees(math.acos(c))
            out[round(th, 9)] = nat
    return out


def lopes_angles(nmax_atoms):
    """Lopes dos Santos, Peres, Castro Neto, PRB 86, 155449 (2012), Eq. (2)-(4): coprime m, r > 0,
    cos(theta) = (3 m^2 + 3 m r + r^2/2) / (3 m^2 + 3 m r + r^2); N = 4 (3 m^2 + 3 m r + r^2) atoms when
    gcd(r, 3) = 1 and N = 4 (m^2 + m r + r^2/3) when gcd(r, 3) = 3."""
    out = {}
    for m in range(1, nmax_atoms):
        if 12 * m * m > nmax_atoms * 3:
            break
        for r in range(1, nmax_atoms):
            if math.gcd(m, r) != 1:
                continue
            s = 3 * m * m + 3 * m * r + r * r
            nat = 4 * s if r % 3 else 4 * (m * m + m * r + r * r // 3)
            if nat > nmax_atoms:
                if r * r > nmax_atoms:   # both formulas grow with r
                    break
                continue
            c = (3 * m * m + 3 * m * r + r * r / 2) / s
            th = math.degrees(math.acos(c))
            out[round(th, 9)] = nat
    return out


def tbg_table(nmax_atoms):
    """Distinct twist angles in (0, 60) degrees of primitive TBG cells (4 N(alpha) atoms <= nmax_atoms), from the
    Hermann search, with the atom count and one index set. Angles theta and -theta and theta + 60 k describe the
    same Bravais coincidence; the TBG angle in (0, 60) is |theta| reduced modulo 60 into that range."""
    out = {}
    for m in search_tbg_hermann(nmax_atoms // 4):
        d = describe(m)
        th = abs(d['theta']) % 60.0
        th = min(th, 60.0 - th) if False else th     # keep (0, 60): theta and 60 - theta are different TBG structures
        if th < 1e-9 or abs(th - 60.0) < 1e-9:
            continue
        key = round(th, 9)
        nat = 4 * d['n1']
        if key not in out or out[key][0] > nat:
            out[key] = (nat, m)
    return out


def legacy_genmoire_cells(nmax_atoms, amax=100, nmoire=None):
    """Replay of the search loops of the old GenMoire.f90 for two identical lattices: a = 0..amax, b = 0..a,
    a' = 0..a, b' = 0..a, same norm (|alpha| = |beta|), the angle compared with the target is its acos (unsigned, in
    [0, 180]); nmoire = (nM1, nM2) also applies its filter nM1 <= det(M - M') <= nM2 (default in the file: 1, 45).
    Returns {angle in (0, 60): smallest atom count found}."""
    out = {}
    for a in range(0, amax + 1):
        for b in range(0, a + 1):
            n = norm((a, b))
            if n == 0 or 4 * n > nmax_atoms:
                continue
            for ap in range(0, a + 1):
                for bp in range(0, a + 1):
                    if norm((ap, bp)) != n:
                        continue
                    if nmoire and not (nmoire[0] <= norm((a - ap, b - bp)) <= nmoire[1]):
                        continue
                    th = legacy_angle((a, b, ap, bp))
                    if th < 1e-9 or th > 60.0 - 1e-9:
                        continue
                    k = round(th, 7)
                    out[k] = min(out.get(k, 10 ** 9), 4 * n)
    return out


# ---------------------------------------------------------------- self-test
def selftest(verbose=True):
    """Checks: the table of the task (angles, atoms, same Bravais cell), the primitivity rule of (a, b, b, a),
    the reduction formula for a - b = 0 mod 3, and Mele = Lopes dos Santos = Hermann for TBG up to 12000 atoms."""
    ok = True

    def check(cond, msg):
        nonlocal ok
        if verbose:
            print(('PASS ' if cond else 'FAIL ') + msg)
        ok = ok and cond

    # The table of the task. Its "primitive indices" are (a+2b, a-b, 2a+b, b-a)/3: same |theta| and atom count, but
    # layer 2 turned the OTHER way (the mirror structure). The same-handedness primitive cell is that formula with the
    # two layer pairs exchanged, ((2a+b)/3, (b-a)/3, (a+2b)/3, (a-b)/3).
    table = [(21.7868, (2, 1, 1, 2), (2, 1, 1, 2), 28, 28),
             (38.2132, (4, 1, 1, 4), (2, 1, 3, -1), 84, 28),
             (27.7958, (5, 2, 2, 5), (3, 1, 4, -1), 156, 52),
             (46.8264, (7, 1, 1, 7), (3, 2, 5, -2), 228, 76),
             (17.8966, (7, 4, 4, 7), (5, 1, 6, -1), 372, 124)]
    for ang, orig, prim, nat_o, nat_p in table:
        do, dp = describe(orig), describe(prim)
        check(abs(abs(do['theta']) - ang) < 5e-5 and abs(abs(dp['theta']) - ang) < 5e-5,
              f'{orig} and {prim}: |theta| {abs(do["theta"]):.4f}, {abs(dp["theta"]):.4f} (table {ang})')
        check(4 * do['n1'] == nat_o and 4 * dp['n1'] == nat_p and dp['primitive'],
              f'{orig}: {4 * do["n1"]} atoms, {prim}: {4 * dp["n1"]} atoms, primitive')
        same = prim if orig == prim else (prim[2], prim[3], prim[0], prim[1])
        if orig != prim:
            check(abs(do['theta'] + dp['theta']) < 1e-9,
                  f'{prim}: layer 2 at {dp["theta"]:+.4f}, opposite to {orig} ({do["theta"]:+.4f}): the mirror cell')
        ds = describe(same)
        q, r = divmod_e((orig[0], orig[1]), (same[0], same[1]))
        q2, r2 = divmod_e((orig[2], orig[3]), (same[2], same[3]))
        check(abs(ds['theta'] - do['theta']) < 1e-9 and r == (0, 0) and r2 == (0, 0) and q == q2
              and same_cell(ds['prim'], do['prim']),
              f'{same}: same signed angle as {orig}, and {orig} = {q} x {same} on both layers (index {norm(q)})')

    # (a, b, b, a), coprime a > b: primitive iff a - b != 0 mod 3; reduction formula otherwise
    bad = 0
    for a in range(1, 60):
        for b in range(1, a):
            if math.gcd(a, b) != 1:
                continue
            d = describe((a, b, b, a))
            if d['primitive'] != ((a - b) % 3 != 0):
                bad += 1
            if (a - b) % 3 == 0:
                p = ((2 * a + b) // 3, (b - a) // 3, (a + 2 * b) // 3, (a - b) // 3)
                dp = describe(p)
                if not (dp['primitive'] and abs(dp['theta'] - d['theta']) < 1e-9 and d['mult'] == 3
                        and same_cell(p, d['prim'])):
                    bad += 1
    check(bad == 0, '(a,b,b,a) primitive iff a-b != 0 mod 3; else (2a+b, b-a, a+2b, a-b)/3 is the primitive cell, '
                    'same signed angle, 3x smaller (coprime 1 <= b < a < 60)')

    # non-coprime and Eisenstein-only common factors
    d = describe((4, 2, 2, 4))
    check(d['mult'] == 4 and d['prim_n1'] == 7, '(4,2,2,4) = 2 x (2,1,1,2): 4 primitive cells')
    d = describe((5, 6, -1, 10))      # = (2 + tau) x (1, -4, 3, -4) up to a unit: integer gcd 1, Eisenstein gcd norm 7
    check(d['mult'] == 7 and d['prim_n1'] == 13 and math.gcd(5, 6, 1, 10) == 1,
          '(5,6,-1,10): integer gcd 1 but 7 primitive cells (a sqrt7 x sqrt7 supercell of the 13-cell)')
    for bad_in in [(0, 0, 1, 2), (1, 2, 0, 0)]:
        try:
            describe(bad_in)
            check(False, f'{bad_in} refused')
        except ValueError:
            check(True, f'{bad_in} refused (degenerate)')
    # zero indices are fine: (5, 0, 0, 5) = 60 degrees, (1, 0, 1, 0) = no twist
    check(abs(describe((5, 0, 0, 5))['theta'] + 60) < 1e-9 and describe((5, 0, 0, 5))['mult'] == 25,
          '(5,0,0,5): layer 2 at -60 deg, 25 copies of the 1x1 cell')

    # TBG completeness against Mele and Lopes dos Santos
    nmax = 12000
    h = tbg_table(nmax)
    me = mele_angles(nmax)
    lo = lopes_angles(nmax)
    ang = lambda d: sorted(d)
    check(len(me) == len(lo), f'Mele {len(me)} angles, Lopes dos Santos {len(lo)} angles in (0, 60), <= {nmax} atoms')
    match = all(abs(x - y) < 1e-7 for x, y in zip(ang(me), ang(lo))) and all(me[k] == lo[k] for k in me)
    check(match, 'Mele and Lopes dos Santos: same angles and atom counts')
    hk = {round(k, 7): v[0] for k, v in h.items()}
    mk = {round(k, 7): v for k, v in me.items()}
    check(set(hk) == set(mk), f'Hermann (signed indices) {len(hk)} angles = Mele {len(mk)} (two-way)')
    check(all(hk[k] == mk[k] for k in mk if k in hk), 'Hermann primitive atom counts = Mele, every angle')
    for filt in (None, (1, 45)):
        old = legacy_genmoire_cells(nmax, nmoire=filt)
        n_found = sum(1 for k in hk if k in old)
        n_prim = sum(1 for k in hk if k in old and old[k] == hk[k])
        filt_txt = ", det(M - M') in [1, 45]" if filt else ''
        if verbose:
            print(f'     old GenMoire loops (indices >= 0, b <= a, a\', b\' <= a, a <= 100'
                  f'{filt_txt}): {n_found} of {len(hk)} angles found with <= {nmax} '
                  f'atoms, {n_prim} with the primitive cell, {len(hk) - n_found} missed')
    return ok


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('indices', nargs='*', type=int, help='m1 m2 m3 m4: describe this cell')
    p.add_argument('--selftest', action='store_true')
    p.add_argument('--tbg-table', type=int, metavar='NATOMS', help='list the primitive TBG cells up to NATOMS atoms')
    a = p.parse_args()
    if a.selftest:
        raise SystemExit(0 if selftest() else 1)
    if a.tbg_table:
        for k, (nat, m) in sorted(tbg_table(a.tbg_table).items()):
            print(f'{k:14.9f} {nat:8d}  {m[0]:5d} {m[1]:5d} {m[2]:5d} {m[3]:5d}')
        return
    if len(a.indices) == 4:
        d = describe(a.indices)
        for k, v in d.items():
            print(f'{k:10s} {v}')


if __name__ == '__main__':
    main()
