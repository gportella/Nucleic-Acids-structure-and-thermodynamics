#include <gromacs/macros.h>
#include <gromacs/smalloc.h>
#include <gromacs/tpxio.h>
#include <gromacs/gmx_fatal.h>
#include <gromacs/typedefs.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define WHAM_MAXFILELEN 2048

typedef struct{
	char *fn;
	double *x[2];
	double *xp[2];
	double *fx[2];
	int    nts;
	double c;
	double tm2;
	double dh2;
	double bf;
	double ef;
	double conc;
} f_info;

void readfiles( const char *fn, f_info *fr, int *nfilesRet); 
