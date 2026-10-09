# TAPW debugging tools

These scripts visualize reciprocal-space debug tables written by GRABNES.
They are development aids, not part of the solver interface.

Install their Python dependencies with:

```sh
python3 -m pip install -r tools/requirements.txt
```

Plot selected Brillouin-zone files:

```sh
python3 tools/tapw/plot_brillouin_zone.py path/to/brillouin_zones_debug_K1.dat \
  path/to/brillouin_zones_debug_K2.dat --g-vectors path/to/g_vectors_debug.dat
```

Summarize all recognized TAPW debug output in a calculation directory:

```sh
python3 tools/tapw/visualize_tapw_debug.py path/to/calculation --output tapw_debug.png
```

Both programs save an image without requiring an interactive display. Add
`--show` when working in a desktop environment. Input and output paths are
explicit, so the tools do not depend on the repository's location.
