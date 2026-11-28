---
short_title: Guide
---

# Welcome to the GRABNES documentation

This documentation covers everything needed to run GRABNES, from installing
the code to interpreting production results.

:::::{grid} 1 2 2 4

::::{grid-item-card} Installation
:link: installation.md

Compile the code, or run it in the browser via Colab.
::::

::::{grid-item-card} Theory
:link: theory/index.md

The physics: tight-binding models, moiré supercells and observables.
::::

::::{grid-item-card} Examples
:link: examples/index.md

The supported, reproducible example calculations.
::::

::::{grid-item-card} Tutorials
:link: tutorials/index.md

The workflows, step by step.
::::

::::{grid-item-card} User guide
:link: user-guide/index.md

The input file and every input key.
::::

::::{grid-item-card} Benchmarks
:link: benchmarks/index.md

Published results.
::::

::::{grid-item-card} Development guide
:link: development/index.md

The state of the code.
::::

::::{grid-item-card} API reference
:link: reference/index.md

The Fortran units.
::::

:::::

## Getting started

GRABNES is a Fortran code. This can be daunting when you have never used this
language, so to give you the opportunity to get a first taste of the GRABNES
code without worrying about compilation issues and compiler questions, we
have compiled a stable version of the code on Google Colab where you can
interact with it using a Jupyter Notebook with Python as the language of
choice. For heavier use of the capabilities of the code, including its OMP
parallelization capabilities, we recommend compiling the code yourself on your
cluster — see [](installation.md).

{button}`Open the interactive tutorial in Google Colab <https://colab.research.google.com/drive/12_ofrjJi4pr6U9j-ER1WocYJkqxEffs-?usp=sharing>`

```{note}
The Colab notebook was set up for the Twistronics 2023 workshop and uses the
key names from before October 2026; its shared resources may no longer be
accessible. To run the current code, build it ([](installation.md)) and start
with the [](examples/index.md).
```

