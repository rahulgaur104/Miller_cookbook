# Miller cookbook
Local equilibrium to GX/GS2 interface

This script creates a local Miller equilibrium for Tokamaks and generates all the relevant data needed for a gyrofluid calculation with GX or a gyrokinetic calculation with GS2 of that equilibrium.

## Limitations
For the moment:
* The code only works with Miller(up-down symmetric) equilibria
* There is no algorithm to handle more than 2 magnetic wells. Non-trivial to add arbitrary wells. See Jessica Baumgaertel's [thesis](https://dataspace.princeton.edu/handle/88435/dsp010r9673776)(FIGG code).
* Only kx = 0 data can be generated. Not hard to generalize.
* No guarantee that the local equilibrium is MHD stable. Difficult.

## Relevant papers
* [Non-circular, finite aspect ratio, local equilibrium model](https://aip.scitation.org/doi/10.1063/1.872666) R. L. Miller et al.
* [Construction of local axisymmetric MHD equilibria](https://inis.iaea.org/search/searchsinglerecord.aspx?recordsFor=SingleRecord&RN=17000660) C.M. Bishop
* [Mercier C and Luc N 1974 Technical Report Commission of the European Communities Report No EUR-5127e 140 Brussels]
