#!/usr/bin/env python3
"""Draft of Section.Name forms for the input keys that have no section (see coding-conventions.md).

Reads docs/user-guide/input-keys.md and writes docs/development/input-key-renaming-proposal.md: one row per
key with the proposed name, grouped by section. The rules below are a first guess from the name of the key
and the routine that reads it; the table is meant to be edited by hand before anything is implemented.

    python3 tools/input/propose_key_names.py
"""
import os
import re

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
SRC = os.path.join(ROOT, "docs", "user-guide", "input-keys.md")
OUT = os.path.join(ROOT, "docs", "development", "input-key-renaming-proposal.md")

# (section, pattern on the lower-case key, number of leading characters of the key dropped from the new name)
RULES = [
    ("Run", r"^(prefix|description|typeofsystem)$", 0),
    ("Output", r"^(writedatafiles|writepos|readdatafiles|keepwavefunction|printbubble)$", 0),
    ("Structure", r"^(xyzfile|sublatticefile|usesublatticefile|readrigidxyz|rigidpositions|readinterlayerdistances|"
                  r"interlayerdistancefile|displacementsfile|invertdisplacements|cellsize|supercell.*|moirecellparameters|"
                  r"interlayerdistance|cellheight|latticeparameter.*|ribbontype|grapheneunitcell|atomsorderdeactivated|"
                  r"basedonmoirecellparameters|singlelayerxyz|bernalreadxyz|readlayerindex|.*width|.*shiftfactor|"
                  r"twistedbladdshift|bernalshift|bridgeshift|trilayeraddshift|createbldomainboundary|bldomainboundary.*)$", 0),
    ("Stack", r"^(onelayer|twolayers.*|threelayers.*|fourlayers.*|fivelayers.*|sixlayers.*|sevenlayers.*|eightlayers.*|"
              r"tenlayers.*|twentylayers.*|encapsulated.*|t2gbn|t2bg|t3bg|t3gwithbn|bnt2gbn|helicaltwistedmbm|"
              r"middletwist|threelayershort|gbntwolayers.*|bnbntwolayers.*|bulk|bulksmall|nonbulksmall)$", 0),
    ("Moire", r"^moire", 5),
    ("TBG", r"^tbg", 3),
    ("GBN", r"^gbn", 3),
    ("BNBN", r"^bnbn", 4),
    ("Intralayer", r"^(f2g2model|koshinointralayer|mayouintralayer|singlelayer.*|bilayert.*|forcebilayerf2g2intralayer|"
                   r"useoldgraphenef2g2|removef2g2flag|deactivateintra.*|deactivateinter?sublattice|deactivate[ab]sublattice)$", 0),
    ("Interlayer", r"^(typeofbl|typeofsl|vpppi0|vppsigma0|bldelta|bilayer.*|coupling.*|renormalize.*|deactivateinterlayer.*|"
                   r"koshinosr|corrugatedinterlayertwocenter|onlyvppsigma|usebng.*|addsecondlayerinteractions|"
                   r"differentcouplings|setnlayerstozero|addpressuredependence|.*parameterset|usetheta.*|findthetasgeometrically|"
                   r"sublattice(in)?dependent|useonlyvab|onlyv0|switchv3sign|deactivatev[36]|newfittingfunctions|"
                   r"addexponentialdecayfordihedral|oppositedxdy|changelatticeparameterforsrivanimodel|interface.*|"
                   r"c[ab][ab]0?|c[ab]p[ab]p|phi[ab][ab]|phi[ab]p[ab]p)$", 0),
    ("Strain", r"^(realstrain.*|onlyfirstneighborrealstrain|periodicstrain|shells.*|strainedmoire|randomstrain|"
               r"realisticbubbles|bigbubble|manybubbles|bubble.*|gaussheight.*)$", 0),
    ("Disorder", r"^(anderson.*|gaussdisorder|gauss.*|delta.*|sublatt(amp|pct)|sublatticedisorder|bubbles|onsiteshift|"
                 r"checker.*)$", 0),
    ("Potential", r"^(sinus.*|cosinus.*|squarefunction.*|squarechecker.*|twodimension.*|armchairshape|addzterm|zterm.*|"
                  r"pnp.*|cdw.*|helicaltwistedmbm_cdw|addsublatticemassterm|sublatticemassterm|onlybottomlayermassterm|"
                  r"addonsiteenergyshift|onsiteenergyshift|uselayerspecificonsiteenergyterms|fourlayeronsiteshifts|"
                  r"fourlayershift\d)$", 0),
    ("SOC", r"^(.*socterm|lambda.*|soclayer.*)$", 0),
    ("Zeeman", r"^(zeemanterm|pseudozeemanterm|spin|spinpolarized|gzeeman.*)$", 0),
    ("Haldane", r"^(haldane.*|paperorientation)$", 0),
    ("Hubbard", r"^(enablescf|scf.*|hubbard.*|aforder.*|forder.*)$", 0),
    ("Kubo", r"^(recursionnumber|numberofenergypoints|epsilon|energymin|energymax|setseed|seedvalue|pdos.*|"
             r"calculateinterval|.*timestep.*|numberoftimesteps)$", 0),
    ("TAPW", r"^(usetapw|usedensematrixtapw)$", 0),
    ("Diag", r"^(sparsediagsolver)$", 0),
    ("MagField", r"^(frankmagneticfield)$", 0),
]


def propose(key):
    low = key.lower()
    for section, pattern, drop in RULES:
        if re.match(pattern, low):
            rest = key[drop:] if drop else key
            return section, f"{section}.{rest[0].upper()}{rest[1:]}"
    return "(to decide)", ""


def main():
    keys = {}
    for line in open(SRC):
        m = re.match(r"\| `([^`]+)` \| ([^|]*) \| `([^`]*)` \| (.*) \|", line)
        if m and "." not in m.group(1):
            keys.setdefault(m.group(1), m.group(4).split(";")[0].strip())
    by = {}
    for key, where in keys.items():
        section, new = propose(key)
        by.setdefault(section, []).append((key, new, where))
    out = ["# Proposed `Section.Name` forms for the keys without a section", "",
           "Draft generated by `tools/input/propose_key_names.py`; edit this table, then implement it with the alias",
           "mechanism described in `coding-conventions.md` (the former names keep working).", "",
           f"{len(keys)} keys. Sections and counts: " + ", ".join(f"{s} {len(v)}" for s, v in sorted(by.items(), key=lambda x: -len(x[1]))) + ".", ""]
    for section in sorted(by, key=lambda s: (s.startswith("("), s)):
        out += [f"## {section}", "", "| Present key | Proposed key | Read in |", "| --- | --- | --- |"]
        out += [f"| `{k}` | {'`' + n + '`' if n else ''} | {w} |" for k, n, w in sorted(by[section], key=lambda r: r[0].lower())]
        out.append("")
    open(OUT, "w").write("\n".join(out))
    print(f"{len(keys)} keys; " + ", ".join(f"{s} {len(v)}" for s, v in sorted(by.items(), key=lambda x: -len(x[1]))))


if __name__ == "__main__":
    main()
