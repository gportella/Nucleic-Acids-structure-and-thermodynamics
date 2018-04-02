
## Determine thermodynamic properties of duplex and triplex DNA from UV-abs
melting curves. 

The thermodynamic properties can be derived from either
multiple curves at different concentration (which is presumably more accurate)
or from a single concentration absorvance curve.

For duplex DNA we assume a two stated equilibrium such as 

D2 <-> D+D

and for a triplex 

TD2 <-> T + D2 <-> T+D+D

The equations coded are for this molecularity, but it is straightforward to
change it to account for any molecularity. Each concentration in the
equilibrium is expressed as a molar fraction and we use the derivative of the
absorbance with respect to the temperature to fit the van't Hoff enthalpy. The
fit is performed as a non-linear least-square fit. 

If the concentrations are available, we propose the following method

1. Each curve is fitted using a DH value that is common for all
2. A second equation, relationg the concentration present in each curve, its
   melting temperature and the  DH, is introduced as an additional linear 
   relation to be fulfilled. 

All equations are adjusted using a non-linear fit.

## How to compile 

A Makefile for the duplex and triplex versions of the programs
is supplied.

The current version requires a compiled non-MPI version of Gromacs 4.6. 
The program uses Gromacs functionality for parsing the arguments, dealing 
with memory and I/O. 
Before compiling, make sure to source this version of Gromacs, e.g. ```source ~/GMX46/bin/GMXRC```.

The program also requires the GSL and GSL-C-Blas libraries to be installed
and available for linking. 

Guillem Portella, Cambridge 15/04/2014

