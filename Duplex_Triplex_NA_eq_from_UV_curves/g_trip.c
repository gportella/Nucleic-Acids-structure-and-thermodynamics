// GPC 10/03/2015 g_tripc v 0.2 
//
//  Re-write python code for computing triplex melting temperatures
// 

#include <gromacs/copyrite.h>
#include <gromacs/macros.h>
#include <gromacs/pbc.h>
#include <gromacs/smalloc.h>
#include <gromacs/statutil.h>
#include <gromacs/gmx_fatal.h>
#include <gromacs/xtcio.h>
#include <gromacs/tpxio.h>
#include <gromacs/mshift.h>
#include <gromacs/vec.h>
#include <gromacs/typedefs.h>
#include <gromacs/gmx_statistics.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <stdlib.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_blas.h>
#include "sg.h"

#define R (0.0019858)

typedef struct 
{
	double tm2; 
	double dh2;
	double tm3;
	double dh3;
	double c;
	double dg3;
} f_p; 

struct fit_params
{
	double k2, k3;
};

struct fit_data
{
	size_t n;
	double *x;
	double *y;
	double tm2;
	double dh2;
};

void do_solve_eq(double *x, int n, f_p p, double *fd, double *ft);
	
void create_fit_curve(double *x[],double *xf[],int n,f_p s)
{
	int i; 
	double *fd; // fraction duplex
	double *ft; // fraction triplex
	snew(fd,n);
	snew(ft,n);
	do_solve_eq(x[0],n,s,fd,ft);
	for (i = 0; i < n; i++)
	{
		//  ( 1 - a2 + c * ( 1- a3) )/(1+c)
		xf[1][i] = ( 1 - fd[i] + s.c * ( 1- ft[i]) )/(1+s.c); 
		xf[0][i] = x[0][i]; 
	}

}

double equilibrium(double a3, void *p)
{
	double ct_2=1;
	struct fit_params *pp = (struct fit_params *) p ;
	return ( ((( pp->k2*pp->k3*ct_2*(1-a3) ) / ( pp->k2*pp->k3*ct_2*(1-a3) + pp->k2 + 1 ) ) - a3)) ; 
}


void do_solve_eq(double *x, int n, f_p p, double *fd, double *ft)
{
	double *k2;
	double *k3;
	int	max_iter=100;
	snew(k2,n);
	snew(k3,n);
	const gsl_root_fsolver_type *solver_type;
	gsl_root_fsolver *solver;
	gsl_function F;
	struct fit_params pf = {0.0, 0.0};

	for(i=0; i<n; i++){ k2[i] = exp((p.dh2/R)* ((1.0/p.tm2) - (1.0/x[i]))) ; }
	for(i=0; i<n; i++){ k3[i] = exp((p.dh3/R)* ((1.0/p.tm3) - (1.0/x[i]))) ; }

	F.function = &equilibrium;

	/* Allocate a bisection solver and set it to use F */
	solver_type = gsl_root_fsolver_brent;
	solver = gsl_root_fsolver_alloc(solver_type);
	double x_lo=0;
	double x_hi=1.0;
	#pragma omp parallel for private(i,solver) firstprivate(F,pf,k2,k3, x_lo, x_hi) 
	for(i=0; i<n; i++){
		int status;
		solver = gsl_root_fsolver_alloc(solver_type);
		double a3=0;
		double a2=0;
		pf.k2 = k2[i];
		pf.k3 = k3[i];
		F.params = &pf;
		x_lo =0.0;
		x_hi =1.0;
		gsl_root_fsolver_set(solver, &F, x_lo, x_hi);
		int iter=0;
		do
		{
			iter++;
			status = gsl_root_fsolver_iterate (solver);
			a3 = gsl_root_fsolver_root (solver);
			x_lo = gsl_root_fsolver_x_lower (solver);
			x_hi = gsl_root_fsolver_x_upper (solver);
			status = gsl_root_test_interval (x_lo, x_hi,0, 0.001);
		//	if (status == GSL_SUCCESS)
		//	printf ("SS %d %f %f %f\n",i, a3, pf.k2,pf.k3);
		//	printf ("%5d [%.7f, %.7f] %.7f %.7f\n", iter, x_lo, x_hi, a2, x_hi-x_lo);
		}
		while (status == GSL_CONTINUE && iter < max_iter);
		a2=( pf.k2*pf.k3*(1-a3) + pf.k2  ) / ( pf.k2*pf.k3*(1-a3) + pf.k2 + 1 ) ;
		fd[i]=a3;
		ft[i]=a2;

	}
	gsl_root_fsolver_free (solver);

	sfree(k2);
	sfree(k3);
}

int eq_fit (const gsl_vector * x, void *data, 
        gsl_vector * f)
{
	// fitting parameters
	double tm3 = gsl_vector_get (x,0);
	double dh3 = gsl_vector_get (x,1);
	double c   = gsl_vector_get (x,2);

	// constants and data to fit
	size_t n = ((struct fit_data *)data)->n;
	double tm2 = ((struct fit_data *)data)->tm2;
	double dh2 = ((struct fit_data *)data)->dh2;
	double *xx = ((struct fit_data *)data)->x;
	double *yy = ((struct fit_data *)data)->y;
	size_t i;
	f_p fpt;

	fpt.tm2=tm2;
	fpt.dh2=dh2;
	fpt.tm3=tm3;
	fpt.dh3=dh3;
	double *fd; // fraction duplex
	double *ft; // fraction triplex
	snew(fd,n);
	snew(ft,n);

	// numerically solving eq 3 from R&C PNAS 1991
	do_solve_eq(xx,n,fpt,fd,ft);

	for (i = 0; i < n; i++)
	{
		//  acc. to R&C PNAS 1991 it should have a diff 
		//  form, but this one from a triplex paper fits better
		//  ( 1 - a2 + c * ( 1- a3) )/(1+c)
		if ( c < 0 ) {c=100000;}
		double Yi = ( 1 - fd[i] + c * ( 1- ft[i]) )/(1+c); 
		gsl_vector_set (f, i, (Yi - yy[i]) );
	}

	return GSL_SUCCESS;

	sfree(fd);
	sfree(ft);
}

f_p do_fit_triplex(double *xx[],int n,f_p p, bool bVerbose)
{
	int i; 
	int iter=0;
	int status;
	f_p fit;
	const size_t pp = 3;
	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;

	gsl_matrix *covar = gsl_matrix_alloc (pp, pp);
	struct fit_data d = { n, xx[0], xx[1], p.tm2, p.dh2};
	gsl_multifit_function_fdf f;
	double x_init[3] = { p.tm3, p.dh3, p.c };
	gsl_vector_view x = gsl_vector_view_array (x_init, pp);
	T = gsl_multifit_fdfsolver_lmsder;
	s = gsl_multifit_fdfsolver_alloc (T, n, pp);

	// copy data to function
	f.f = &eq_fit;
	f.df = NULL;
	f.fdf = NULL;
	f.p = pp;
	f.n = n;
	f.params = &d;  

	gsl_multifit_fdfsolver_set (s, &f, &x.vector);

	do
	{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);
		if(bVerbose){
			printf ("iter: %3u x = % 15.8f % 15.8f % 15.8f "
					"|f(x)| = %g\n",iter,
					gsl_vector_get (s->x, 0), 
					gsl_vector_get (s->x, 1),
					gsl_vector_get (s->x, 2), 
					gsl_blas_dnrm2 (s->f));
		}

		if (status)
			break;
		status = gsl_multifit_test_delta (s->dx, s->x,
				1e-8, 1e-8);
	}
	while (status == GSL_CONTINUE && iter < 1000);
	gsl_multifit_covar (s->J, 0.0, covar);

	#define FIT(i) gsl_vector_get(s->x, i)
	#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

	fit.tm2 = p.tm2;
	fit.dh2 = p.dh2;
	fit.tm3 = gsl_vector_get(s->x, 0);
	fit.dh3 = gsl_vector_get(s->x, 1);
	fit.c   = gsl_vector_get(s->x, 2);

		double chi = gsl_blas_dnrm2(s->f);
		double dof = n - pp;
		double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 
	if(bVerbose)
	{ 

		printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);
		printf ("Tm3    = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
		printf ("DH3    = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
		printf ("c      = %.5f +/- %.5f\n", FIT(2), c*ERR(2));
		printf ("status = %s\n", gsl_strerror (status));
	}else{

		printf ("Tm3    = %.5f  %.5f\n", FIT(0)-273.15, c*ERR(0));
	}

	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);

		return fit;
}

void scale_data(double *x[],int n, double t_m)
{
	int i;
	double min=0;
	min = x[1][0];
	for(i=0; i<n; i++){
		x[0][i] += 273.15;
		x[1][i] = (x[1][i] - min) / (t_m - min)   ;
	}

}

double find_t_height(double *x[],int n,double m)
{
	int i; 
	double val=0;
	val=x[1][n-1];
	for(i=0; i<n; i++){
		if ( x[1][i] <= val ) { val = x[1][i] ; } 
	}
	val = m - val ;
	return val;
}

double find_maxabs(double *x[],int n,double t)
{
	int i;
	double val=0;
	for(i=0; i<n; i++){
		if ( x[0][i] <= t) { val = x[1][i] ; } 
	}
	return val; 
}

void compute_verivatives(double *x[], double *xp[], int nts)
{
	int i=0;
	double x0=0;
	double y0=0;

	x0=x[0][0];
	y0=x[1][0];
	for (i=1; i<nts; i++){
		xp[0][i-1] = x0+(x[0][i] -x0)/2 ; 
		xp[1][i-1] = (x[1][i] -y0)/(x[0][i]-x0) ;
		x0=x[0][i] ; 
		y0=x[1][i] ;
	}
}


double suggest_tm3(double *ctrip[],int nts, double dT)
{
	int i,j,k;
	double tm3=0;
	double val=0;
	double *dctrip[2];
	double *smth[2];
	// left/right number of frames to fullfill a local maxima
	int PAD=10;
	int lmax=0;
	int rmax=0;
	bool *bMax;
	snew(dctrip[0],nts-1);
	snew(dctrip[1],nts-1);
	snew(smth[0],nts-1);
	snew(smth[1],nts-1);
	snew(bMax,nts-1);
	// idk if they have been initialized to 0, but just in case
	for(i=0; i< nts-1; i++){ bMax[i]=FALSE;}

	compute_verivatives(ctrip,dctrip,nts);
	// smooth a few times
	sgsmooth(dctrip,nts,smth); 
	sgsmooth(smth,nts,smth); 
	sgsmooth(smth,nts,smth); 
	sgsmooth(smth,nts,smth); 
	// search for local maxima within a window of 2*PAD
	// we addjust the padding based on the number of points
	PAD = (int) (nts / 10) ; 
	#pragma omp parallel for private(j,lmax,rmax)
	for(i=PAD; i< nts-PAD-1; i++){
		for(j=1; j<PAD; j++){
			if (smth[1][i-j] < smth[1][i] ) { lmax++; }
			if (smth[1][i+j] < smth[1][i] ) { rmax++; } 
		}
		if(lmax+rmax==2*(PAD-1)) { bMax[i] = TRUE ; }
		lmax=0;
		rmax=0;
	}
	//for(i=0; i< nts-1; i++){ printf("%f %f %f \n",smth[0][i],(float)0.001*bMax[i], dT);}
	// throw away after the duplex melting curve, keep whatever is first
	tm3=0;
	for(i=PAD; i< nts-PAD-1; i++){ 
		if (smth[0][i] < dT && bMax[i] && smth[1][i] > val) {
			val=smth[1][i];
			tm3=smth[0][i];
		}
	}

	return tm3;
}


int main(int argc, char *argv[])
{
	const char         *desc[] = {
		"Read triplex melting curve and duplex data to fit equilibrium [BR]"
			"Solving the equations takes time."
	};

	FILE	*f_melt;
	FILE	*f_fit;
	static real		dT=70;
	static real		dH=-60;
	static real		dth=0.001;
	static real		bf=25;
	static real		ef=100;
	static real	tmaxabs=-1;
	static bool	bVerbose=FALSE;

	/* Command-line arguments */
	t_pargs          pa[] = {
		{"-dm", FALSE, etREAL, {&dT},
			"Duplex melting temperature"},
		{"-dh", FALSE, etREAL, {&dH},
			"Duplex enthalpy"},
		{"-dth", FALSE, etREAL, {&dth},
			"Duplex transition height"},
		{"-bf", FALSE, etREAL, {&bf},
			"Begin fitting triplex curve"},
		{"-ef", FALSE, etREAL, {&ef},
			"End fitting triplex curve"},
		{"-tmax", FALSE, etREAL, {&tmaxabs},
			"T of max absorbance, use end fitting if not given"},
		{"-v", TRUE, etBOOL, {&bVerbose},
			"Be verbose."},
	};

	/* Output files */
	t_filenm            fnm[] = {
		{ efDAT, "-d", "triplex_melt.dat", ffREAD },
		{ efXVG, "-o", "fitted_melt", ffWRITE  },
	};

#define NFILE asize(fnm)
    output_env_t    oenv;
	int				i,j,k;
	int				npargs;
	int				n=0;
	int				nts=0;
	char			buf[256];
	const char		*tripfile=NULL;
	const char		*fitfile=NULL;
	const char		*out_file=NULL; 
	char			line[255];
	double			*ctrip[2]; 
	double			*ftrip[2];
	double			*dc[2]; 
	double			xx,yy;
	double			tm3;
	double			max_abs=0;
	double			t_height=0;
	double			c_fix=0;
	double			dh3=0;
	bool			bsettmax;
	f_p				p;
	f_p				solution; 

	npargs=asize(pa);
//	CopyRight(stderr, argv[0]);
	parse_common_args(&argc, argv, 0,
			NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);

	p.dh2= (double) dH;
	p.tm2= (double) dT+273.15;
	// if we don't set a tmax, use the same as ef
	if(tmaxabs<=0) { tmaxabs = ef ; }
	tripfile = opt2fn("-d",NFILE,fnm);
	if(tripfile) {f_melt = ffopen(tripfile,"r");}
	while(fgets(line, 255, f_melt) != NULL) { 
		sscanf(line, "%lf %lf", &xx, &yy);
		if( xx > bf && xx < ef) {nts++; }
	}
	snew(ctrip[0],nts);
	snew(ctrip[1],nts);
	snew(ftrip[0],nts);
	snew(ftrip[1],nts);
	if (bVerbose) {
		printf("There are %d T from %s\n",nts,tripfile);
	}
	rewind(f_melt);
	// Crude parsing, no comments or extra lines allowed
	j=0;
	xx=0; yy=0;
	while(fgets(line, 255, f_melt) != NULL)
	{
		sscanf(line, "%lf %lf", &xx, &yy);
		if( xx > bf && xx < ef) {ctrip[0][j]=xx; ctrip[1][j]=yy; j++; }
	}
	max_abs=find_maxabs(ctrip,nts,tmaxabs);
	t_height=find_t_height(ctrip,nts,max_abs);
	c_fix = (dth)/ (t_height - dth ) ;
	if (c_fix<0) { 
		c_fix=1; 
		printf("WARNING! Wrong with initial transition height\n");
		printf("WARNING! Setting the ratio of transitions to 1 to proceed.\n");
	} 
	scale_data(ctrip,nts,max_abs);
	p.tm3=suggest_tm3(ctrip,nts,p.tm2);
	p.dh3=-80;
	p.c = c_fix;
	if(bVerbose){
		printf("Abs. max is %5.3f\n",max_abs); 
		printf("We suggest %5.3f as the triplex melting point\n",p.tm3); 
		printf("Transition height is %5.3f, ratio %5.3f \n",t_height, c_fix); 
	}
	//do the fit
	solution=do_fit_triplex(ctrip,nts,p,bVerbose);
	create_fit_curve(ctrip,ftrip,nts,solution);

	// Output experimental vs fitted
	fitfile = opt2fn("-o",NFILE,fnm);
	if(fitfile) {
		f_fit = ffopen(fitfile,"w");
		for(i=0; i<nts; i++){  
			fprintf(f_fit,"%f %f \n", ctrip[0][i]-273.15, ctrip[1][i]);
		}
			fprintf(f_fit,"&\n");
		for(i=0; i<nts; i++){  
			fprintf(f_fit,"%f %f \n", ftrip[0][i]-273.15, ftrip[1][i]);
		}
		fclose(f_fit);
	}

	return 0;

}

