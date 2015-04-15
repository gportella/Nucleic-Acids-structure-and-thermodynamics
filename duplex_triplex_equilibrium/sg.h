#include <stdio.h>
#include <gromacs/typedefs.h>
#include <gromacs/smalloc.h>
#include <math.h>

#define  NP      50    //Maximum number of filter coefficients
#define  MMAX     6    //Maximum order of smoothing polynomial
#define NMX  100

//Note: ind 0 not used here.
typedef  float MAT[MMAX+2][MMAX+2];

float *smoothed, *ysave;
float c[NP+1];
int   ind[51];

int   i,j,m,ndata,nl,nr;
float dt,t,tbegin,temp,tend;

FILE  *fp_in, *fp_out;

int IMin(int ia, int ib);
void LUDCMP(MAT A, int N, int np, int *INDX, int *D, int *CODE);
void LUBKSB(MAT A, int N, int np, int *INDX, float *B);
void savgol(float *c, int np, int nl, int nr, int ld, int m);
void sgsmooth(double *x[],int nts, double *smth[]);

