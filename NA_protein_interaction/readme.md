
## Compute interactions between a group of atoms against two other groups of atoms from MD trajectories (Gromacs XTC files)

Tool to compute both electrostatic and Lennard-Jones interactions, it can
use a cut-off and neigbour searching to speed up the calculation. 
Compare to using gmx mdrun -rerun is slower, but it is a bit more flexible
in defining the interaction groups. 


## How to compile

> Make sure you have sourced Gromacs before you attempt to compile and/or run the program. 
It was originally compiled against Gromacs 2016, but maybe it works with newer versions as well.

To compile, run `cmake` in the directory witht the sources. You might need to help `cmake` a 
bit by providing a path to your C and C++ compilers. You should use the same compilers and libraries
as the ones you used to compile Gromacs.

```bash
cmake . -DCMAKE_C_COMPILER=/opt/local/bin/clang-mp-7.0 -DCMAKE_CXX_COMPILER=/opt/local/bin/clang++-mp-7.0 -DCMAKE_PREFIX_PATH=/opt/local/lib/
make -j 4
```

It compiles as standalone application, it is not merged into the gmx "name_tool". 

Once compiled ./interact -h should give you the options. These should be documented 
fine, but just in case you can also check the code. 





