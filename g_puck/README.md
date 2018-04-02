
## Computes puckering angle and amplitude for ribose DNA/RNA rings from a MD trajectory (Gromacs XTC)

Should be linked against Gromacs-4.6.x (compiled without MPI support).
For convenience source gromacs, e.g. ```source ~/GMX/4.6/bin/GMXRC```, but one could
also edit the Makefile.

The tools computes puckering angle and amplitude for ribose DNA/RNA rings from a Gromacs XTC file.
It outpus the (circular) average of the puckering angle. Beside a trajectory and a gromacs 
topology/PDB file, it needs a file specifying the nucleic acid residues that you are interested in.

Guillem Portella  29/04/2015 g_puckering v 0.1
