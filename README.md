# Miller cookbook
Local equilibrium to GX/GS2 interface

This script creates a local Miller equilibrium for Tokamaks and generates all the relevant data needed for a gyrofluid calculation with GX or a gyrokinetic calculation with GS2 of that equilibrium.

## Generating eikcoefs
Set the Miller parameters in the first 50 lines of the script local_eikcoefs_gen_norm.py and run it. The coefficient files are saved in the directory output_grid_files.

## Limitations
For the moment:
* the code only works with Miller(up-down symmetric) equilibria
* there is no algorithm to handle more than 2 magnetic wells. Non-trivial to add arbitrary wells. See Jessica Baumgaertel's [thesis](https://dataspace.princeton.edu/handle/88435/dsp010r9673776)(FIGG code).
* there is no guarantee that the local equilibrium is MHD stable. Difficult.

## Relevant papers
* [Non-circular, finite aspect ratio, local equilibrium model](https://aip.scitation.org/doi/10.1063/1.872666) R. L. Miller et al.
* [Construction of local axisymmetric MHD equilibria](https://inis.iaea.org/search/searchsinglerecord.aspx?recordsFor=SingleRecord&RN=17000660) C.M. Bishop
* [Mercier C and Luc N 1974 Technical Report Commission of the European Communities Report No EUR-5127e 140 Brussels]
