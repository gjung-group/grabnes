Compile this code with your favorite fortran compiler after modifying some of the input values in it (angle, lattice vectors, optimization scheme and tolerance).

Guide:
phideg, is the desired angle in degrees. The program will try to find the best fit to this angle.
degtol, is the allowed tolerance for the angle. Angles within this tolerance are considered as options by the program.
lattol, tolerance related to the lattice parameter of BN (or the second graphene layer).
aG and aBN0, lattice parameters of graphene and BN respectively. The final lattice parameter of BN can be changed to fit both layers.
select, determines the criteria to select a "best" option from the allowed moires. 'delta' means that you want to keep the lattice parameter as close as possible to the one given by aBN0. 'angle' means you prefer to give more weight to the angle.
If fdfbuild is true, the "best" fit will be written into a fdf Siesta file that you can visualize with GDIS.

Note: GDIS does not seem to be maintained anymore.

Other ways to generate twisted bilayer systems:
- Code by Manish Group (called Twister)
- The small piece of code at the beginning of the PyBinding notebook on twisted bilayer systems (the `PyBinding/` notebooks are no longer part of this repository; see its Git history).

The PyBinding notebook can read any of these systems as long as they are in .xyz format where the usually obsolete second line is used to pass the cell definition information.
