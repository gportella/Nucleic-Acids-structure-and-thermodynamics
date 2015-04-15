#include <stdio.h>
#include <math.h>

#define  NMAX  5048    //Maximum number of input data ordinates
#define  NP      50    //Maximum number of filter coefficients
#define  MMAX     6    //Maximum order of smoothing polynomial

//Note: ind 0 not used here.
typedef  float MAT[MMAX+2][MMAX+2];

float signal[NMAX+1], ysave[NMAX+1];
float c[NP+1];
//int   ind[NP+1];
int   ind[51];

int   i,j,m,ndata,nl,nr;
float dt,t,tbegin,temp,tend;

FILE  *fp_in, *fp_out;


int IMin(int ia, int ib) {
  if (ia<=ib) return ia;
  else return ib;
}


  void LUDCMP(MAT A, int N, int np, int *INDX, int *D, int *CODE);
  void LUBKSB(MAT A, int N, int np, int *INDX, float *B);


  void savgol(float *c, int np, int nl, int nr, int ld, int m)  {
/*------------------------------------------------------------------------------------------- 
 USES lubksb,ludcmp given below. 
 Returns in c(np), in wrap-around order (see reference) consistent with the argument respns 
 in routine convlv, a set of Savitzky-Golay filter coefficients. nl is the number of leftward 
 (past) data points used, while nr is the number of rightward (future) data points, making 
 the total number of data points used nl +nr+1. ld is the order of the derivative desired 
 (e.g., ld = 0 for smoothed function). m is the order of the smoothing polynomial, also 
 equal to the highest conserved moment; usual values are m = 2 or m = 4. 
-------------------------------------------------------------------------------------------*/
  int d,icode,imj,ipj,j,k,kk,mm;
  int indx[MMAX+2];
  float fac, sum;
  MAT   a;
  float b[MMAX+2];

  if (np<nl+nr+1 || nl<0 || nr<0 || ld>m || m>MMAX || nl+nr<m) {
    printf("\n Bad args in savgol.\n");
    return;
  }

  for (i=1; i<=MMAX+1; i++) {
    for (j=1; j<=MMAX+1; j++) a[i][j]=0.0;
    b[i]=0.0;
    indx[i]=0;
  }

  for (ipj=0; ipj<=2*m; ipj++) { //Set up the normal equations of the desired leastsquares fit.
    sum=0.0;
    if (ipj==0) sum=1.0;
    for (k=1; k<=nr; k++) sum += (float) pow(k,ipj);
    for (k=1; k<=nl; k++) sum += (float) pow(-k,ipj);
    mm=IMin(ipj,2*m-ipj);
    imj=-mm;
    do {
      a[1+(ipj+imj)/2][1+(ipj-imj)/2]=sum;
      imj += 2;
    } while (imj<=mm); 
  }

  LUDCMP(a,m+1,MMAX+1,indx,&d,&icode);    //Solve them: LU decomposition
 
  for (j=1; j<=m+1; j++) b[j]=0.0;
  b[ld+1]=1.0;    //Right-hand side vector is unit vector, depending on which derivative we want.

  LUBKSB(a,m+1,MMAX+1,indx,b);      //Backsubstitute, giving one row of the inverse matrix.
  
  for (kk=1; kk<=np; kk++)          //Zero the output array (it may be bigger than the number
    c[kk]=0.0;                      //of coefficients.

  for (k=-nl; k<=nr; k++) {         //Each Savitzky-Golay coefficient is the dot product
    sum=b[1];                       //of powers of an integer with the inverse matrix row.
    fac=1.0;
    for (mm=1; mm<=m; mm++) {
      fac *= k;
      sum += b[mm+1]*fac;
    } 
    kk=((np-k) % np) + 1;           //Store in wrap-around order}
    c[kk]=sum;
  }
}

/**************************************************************
* Given an N x N matrix A, this routine replaces it by the LU *
* decomposition of a rowwise permutation of itself. A and N   *
* are input. INDX is an output vector which records the row   *
* permutation effected by the partial pivoting; D is output   *
* as -1 or 1, depending on whether the number of row inter-   *
* changes was even or odd, respectively. This routine is used *
* in combination with LUBKSB to solve linear equations or to  *
* invert a matrix. Return code is 1, if matrix is singular.   *
**************************************************************/
void LUDCMP(MAT A, int N, int np, int *INDX, int *D, int *CODE) {

#define NMX  100

 float AMAX,DUM,SUM,TINY;
 float VV[NMX];
 int   I,IMAX,J,K;

 TINY = (float)1e-12;

 *D=1; *CODE=0;

 for (I=1; I<=N; I++) {
   AMAX=0.0;
   for (J=1; J<=N; J++)
     if (fabs(A[I][J]) > AMAX) AMAX=(float) fabs(A[I][J]);
   if (AMAX < TINY) {
     *CODE = 1;
     return;
   }
   VV[I] = (float)1.0 / AMAX;
 }

 for (J=1; J<=N; J++) {
   for (I=1; I<J; I++) {
     SUM = A[I][J];
     for (K=1; K<I; K++)
       SUM -= A[I][K]*A[K][J];
     A[I][J] = SUM;
   }
   AMAX = 0.0;
   for (I=J; I<=N; I++) {
     SUM = A[I][J];
     for (K=1; K<J; K++)
       SUM -= A[I][K]*A[K][J];
     A[I][J] = SUM;
     DUM = VV[I]*(float)fabs(SUM);
     if (DUM >= AMAX) {
       IMAX = I;
       AMAX = DUM;
     }
   }  
   
   if (J != IMAX) {
     for (K=1; K<=N; K++) {
       DUM = A[IMAX][K];
       A[IMAX][K] = A[J][K];
       A[J][K] = DUM;
     }
     *D = -(*D);
     VV[IMAX] = VV[J];
   }

   INDX[J] = IMAX;
   if ((float)fabs(A[J][J]) < TINY)  A[J][J] = TINY;

   if (J != N) {
     DUM = (float)1.0 / A[J][J];
     for (I=J+1; I<=N; I++)  A[I][J] *= DUM;
   } 
 } //j loop

} //LUDCMP()


/*****************************************************************
* Solves the set of N linear equations A . X = B.  Here A is    *
* input, not as the matrix A but rather as its LU decomposition, *
* determined by the routine LUDCMP. INDX is input as the permuta-*
* tion vector returned by LUDCMP. B is input as the right-hand   *
* side vector B, and returns with the solution vector X. A, N and*
* INDX are not modified by this routine and can be used for suc- *
* cessive calls with different right-hand sides. This routine is *
* also efficient for plain matrix inversion.                     *
*****************************************************************/
void LUBKSB(MAT A, int N, int np, int *INDX, float *B)  {

  float SUM;
  int I,II,J,LL;

  II = 0;

  for (I=1; I<=N; I++) {
    LL = INDX[I];
    SUM = B[LL];
    B[LL] = B[I];
    if (II != 0)
      for (J=II; J<I; J++)
        SUM -= A[I][J]*B[J];
    else if (SUM != 0.0)
      II = I;
    B[I] = SUM;
  }

  for (I=N; I>0; I--) {
    SUM = B[I];
    if (I < N)
      for (J=I+1; J<=N; J++)
        SUM -= A[I][J]*B[J];
    B[I] = SUM / A[I][I];
  }

}


void main()  {

  for (i=1; i<=NMAX; i++) signal[i]=0.0;

  //open input and output file
  fp_in=fopen("smooth.dat","r");
  fp_out=fopen("tsavgol.lst","w");

  //read number of input signal points in input file
  fscanf(fp_in, "%d", &ndata);

  //read ndata couples T(i), Y(i) in input data file
  for (i=1; i<=ndata; i++) {
    fscanf(fp_in, "%f %f", &temp, &signal[i]);
    if (i == 1) tbegin=temp;
    if (i == ndata) tend=temp;
  }
  fclose(fp_in);

  for (i=1; i<=NMAX; i++) ysave[i]=signal[i]; //save unsmoothed signal

  nl=8; nr=8; m=4;                            //see savgol
  
// seek shift ind for given case nl, nr, m (see savgol).
  ind[1]=0;
// example: case nl=nr=5
// ind(2)=-1; ind(3)=-2; ind(4)=-3; ind(5)=-4; ind(6)=-5
  j=3;
  for (i=2; i<=nl+1; i++) {
    ind[i]=i-j;
    j += 2;
  }
// ind(7)= 5; ind(8)= 4; ind(9)= 3; ind(10)=2; ind(11)=1
  j=2;
  for (i=nl+2; i<=nl+nr+1; i++) {
    ind[i]=i-j;
    j += 2;
  }

// calculate Savitzky-Golay filter coefficients.
  savgol(c,nl+nr+1,nl,nr,0,m);

  printf("\n Number of left points .......: %d\n", nl);
  printf(" Number of right points ......: %d\n", nr);
  printf(" Order of smoothing polynomial: %d\n\n", m);
  printf(" Savitzky-Golay Filter Coefficients:\n");
  for (i=1; i<=nl+nr+1; i++) printf("%10.6f", c[i]);
  printf("\n");

// Apply filter to input data.
  for (i=1; i<=ndata-nr; i++) {
    signal[i]=0.0;
    for (j=1; j<=nl+nr+1; j++)
	  if (i+ind[j]>0)   //skip left points that do not exist.
        signal[i] += c[j]*ysave[i+ind[j]];
  }

// write results to output file
  dt=(tend-tbegin)/(ndata-1);
  t=tbegin-dt;
  for (i=1; i<=ndata; i++) {
    t += dt;
    fprintf(fp_out,"   %10.6f   %10.6f   %10.6f\n", t, ysave[i], signal[i]);
  }
  
  fclose(fp_out);
  printf("\n Results in file tsavgol.lst.\n\n");

}
