# Mercier-Luc-local-cookbook
Local equilibrium to GS2 interface

This script creates a local Miller equilibrium for Tokamaks and generates all the relevant data needed for a gyrofluid calculation with GX or gyrokinetic calculation with GS2 of that equilibrium.

## Limitations
For the moment:
* The code only works with Miller equilibria
* There is no algorithm to handle more than 2 magnetic wells
* Only kx = 0 data can be generated 
* No guarantee that the local equilibrium is MHD stable

## Relevant papers:\
* [Non-circular, finite aspect ratio, local equilibrium model](https://aip.scitation.org/doi/10.1063/1.872666)R. L. Miller et al.\\
* [Construction of local axisymmetric MHD equilibria](https://inis.iaea.org/search/searchsinglerecord.aspx?recordsFor=SingleRecord&RN=17000660)C.M. Bishop\\
* [Mercier C and Luc N 1974 Technical Report Commission of the European Communities Report No EUR-5127e 140 Brussels]
