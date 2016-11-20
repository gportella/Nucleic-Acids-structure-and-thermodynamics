
Compute interactions between on group of atoms against
two other groups of atoms from MD trajectories (Gromacs XTC files)
It computes both electrostatic and Lennard-Jones interactions, it can
use a cut-off and neigbour searching to speed up the calculation. 
Compare to using gmx mdrun -rerun is slower, but it is a bid more flexible
in defining the interaction groups. 


Requires Gromacs 2016

It compiles as standalone application, it is not merged into the gmx "name_tool" 
mechanism for convenience. 

Once compiled ./interact -h should give you the options. These should be documented 
fine, but just in case you can also check the code. 





