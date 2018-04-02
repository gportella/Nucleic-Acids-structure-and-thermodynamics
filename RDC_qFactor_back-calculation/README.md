
## Back-calculate RDCs from MD trajectories (Gromacs XTC files) using the theta method. 

Performs a rotational search for the optimal correlation between
calculated and experimental RDCs, and then back-calculates the 
RDCs. The calculate RDCs average both in time and across an ensemble. 


## How to compile and use

**Important**: requires Gromacs 4.6.3 compiled without mpi support, aka serial version.

After sourcing Gromacs 4.6.3, e.g. ```source ~/gmx463/bin/GMXRC``` just issue ```make``` to compile.

Once compiled ./g_qrdc -h outputs the options to run the program. These should be documented 
fine, but just in case you can also check the code. 




