// GPC 10/03/2015 g_tripc v 0.2 
//
//  Re-write python code for computing triplex melting temperatures
//  Each curve is fitted separately, and then we compute the r2 for
//  the linear relationship 1/tm vs ln(Ct). We run N iterations
//  and keep the best one
//
//  TODO Would be nice to clean up the code a bit, specially the part
//  were I do multiple triplex fittings, should be a call to a function
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
#include "read_mult_t_inp.h"

#define R (0.0019858)
struct fit_params
{
    double k2, k3;
};

typedef struct{
	int max_MC_iter;
	real wt;
} fit_option;

struct fit_data
{
    size_t n;
    int  w;
	f_info *f;
};

typedef struct{
	double dh3;
	double ds3;
	double dg3;
	double r2;
} f_all_sol; 

typedef struct
{
	double tm2;
	double dh2;
	double tm3;
	double dh3;
	double c;
	double dg3;
} f_p;

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
    int max_iter=200;
    int status;
    int iter=0;
    snew(k2,n);
    snew(k3,n);
    const gsl_root_fsolver_type *solver_type;
    gsl_root_fsolver *solver;
    gsl_function F;
    struct fit_params pf = {0.0, 0.0};


	// van't hoff for each equilibrium (duplex and triplex)
    for(i=0; i<n; i++){ k2[i] = exp((p.dh2/R)* ((1.0/p.tm2) - (1.0/x[i]))) ; }
    for(i=0; i<n; i++){ k3[i] = exp((p.dh3/R)* ((1.0/p.tm3) - (1.0/x[i]))) ; }

    F.function = &equilibrium;

    /* Allocate a bisection solver and set it to use F */
    solver_type = gsl_root_fsolver_brent;
    solver = gsl_root_fsolver_alloc(solver_type);
    double x_lo=0;
    double x_hi=1.0;
    for(i=0; i<n; i++){
       // solver = gsl_root_fsolver_alloc(solver_type);
        double a3=0;
        double a2=0;
        pf.k2 = k2[i];
        pf.k3 = k3[i];
        F.params = &pf;
        x_lo =0.0;
        x_hi =1.0;
        gsl_root_fsolver_set(solver, &F, x_lo, x_hi);
		iter=0;
        do
        {
            iter++;
            status = gsl_root_fsolver_iterate (solver);
            a3 = gsl_root_fsolver_root (solver);
            x_lo = gsl_root_fsolver_x_lower (solver);
            x_hi = gsl_root_fsolver_x_upper (solver);
            status = gsl_root_test_interval (x_lo, x_hi,0, 0.0001);
          if (status == GSL_SUCCESS){
         // printf ("SS %d %f %f %f\n",i, a3, pf.k2,pf.k3);
         // printf ("%5d [%.7f, %.7f] %.7f %.7f\n", iter, x_lo, x_hi, a2, x_hi-x_lo);
		  }
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

// function that the gls solver calls to make the fit
// for each individual curve
int eq_fit_trip(const gsl_vector * x, void *d,
        gsl_vector *fit)
{
	int i,j;
	size_t n = ((struct fit_data *)d)->n;
	int  w = ((struct fit_data *)d)->w;
	double *fd; // fraction duplex
	double *ft; // fraction triplex
	int ntot=0;
	// fitting parameters
	double tm3 = gsl_vector_get (x,0);
	double c = gsl_vector_get (x,1);
	double dh3   = gsl_vector_get (x,2);
	// constants and data to fit
	double tm2 = ((struct fit_data *)d)->f[w].tm2;
	double dh2 = ((struct fit_data *)d)->f[w].dh2;
	double nts = ((struct fit_data *)d)->f[w].nts;
	double *xx = ((struct fit_data *)d)->f[w].x[0];
	double *yy = ((struct fit_data *)d)->f[w].x[1];

	f_p fpt;
	fpt.tm2=tm2;
	fpt.dh2=dh2;
	fpt.tm3=tm3;
	fpt.dh3=dh3;
	snew(fd,nts);
	snew(ft,nts);

	// numerically solving eq 3 from R&C PNAS 1991
	do_solve_eq(xx,nts,fpt,fd,ft);

	for (j = 0; j < nts; j++)
	{
		//  acc. to R&C PNAS 1991 it should have a diff 
		//  form, but this one from a triplex paper fits better
		//  ( 1 - a2 + c * ( 1- a3) )/(1+c)
		if ( c < 0 ) {c=100000;}
		double Yi = ( 1 - fd[j] + c * ( 1- ft[j]) )/(1+c);
		gsl_vector_set (fit, j, (Yi - yy[j]) );
		ntot++;
	}
	sfree(fd);
	sfree(ft);

	return GSL_SUCCESS;

}

// function that the gls solver calls to make the fit
// for 1/tm3 vs ln(Ct)l 
int eq_fit_straight(const gsl_vector * x, void *d,
        gsl_vector *fit)
{
	int i;
	double dh3 = gsl_vector_get (x,0);
	double ds3 = gsl_vector_get (x,1);
	size_t n = ((struct fit_data *)d)->n;
	for (i = 0; i < n; i++){
		double conc = ((struct fit_data *)d)->f[i].conc;
		double tm3 = ((struct fit_data *)d)->f[i].tm3;
		gsl_vector_set (fit,i, 
				(1.0/tm3  - ((R/dh3)*log(conc) + ((ds3-R*log(4))/dh3)) )  );
	}

	return GSL_SUCCESS;
}

void eval_fit_straight(f_all_sol mf, f_info *f_inp, int nf, real *f)
{
	int i;
	for(j=0; j<nf; j++){
		f[j] = (R/mf.dh3)*log(f_inp[j].conc) + ((mf.ds3-R*log(4))/mf.dh3);
	}

}

f_all_sol do_fit_triplex(f_info *f,int n, real *t_min, real *t_max, bool bVerbose, bool bUpdate)
{
	int i,j,k;
	int iter=0;
	int tot_points=0;
	int status;
	int tot_abs=0;
	int tot_inv_t=0;
	real inv_t=0;
	real tss_inv=0;
	real mean_inv_t=0;
	real abs_t=0;
	real mean_abs_t=0;
	real tss_abs=0;
	real tss_tot=0;
	f_all_sol fit;
	f_info *local_f;
	size_t pp = 0;

	snew(local_f,n);

	// count points within boundaries, allocate 
	for(i=0; i<n; i++){
		for(j=0; j<f[i].nts; j++){
			if( f[i].x[0][j] >= t_min[i] && f[i].x[0][j] <= t_max[i] ) { local_f[i].nts++ ; }
		}
		snew(local_f[i].x[0],local_f[i].nts);
		snew(local_f[i].x[1],local_f[i].nts);
	}
	// now copy the data
	// also the conc
	k=0;
	for(i=0; i<n; i++){
		local_f[i].conc = f[i].conc ;
		local_f[i].dh2 = f[i].dh2 ;
		local_f[i].tm2 = f[i].tm2 ;
		for(j=0; j<f[i].nts; j++){
			if( f[i].x[0][j] >= t_min[i] && f[i].x[0][j] <= t_max[i] ) {
				local_f[i].x[0][k] = f[i].x[0][j] ;
				local_f[i].x[1][k] = f[i].x[1][j] ;
				// we use this loop to compute the TSS, the total sum of squares of the
				// "y" data, to be used for r_squared after we know chi_sq
				abs_t += local_f[i].x[1][k] ;
				tot_abs++;
				k++;
			}
		}
		k=0;
	}

	// the number of parameters to fit is 3
	// corresponding to tm3, dh3 and c, 
	pp = 3;
		const gsl_multifit_fdfsolver_type *T;
		T = gsl_multifit_fdfsolver_lmsder;
		gsl_multifit_fdfsolver *s;
	// do a fit for each triplex curve
	for (i=0; i<n; i++){
		// Total number of points to fit
		if (bVerbose){printf("Working on curve n %d\n",i);} 
		struct fit_data d = { n,i, local_f};
		tot_points = local_f[i].nts ;
		gsl_matrix *covar = gsl_matrix_alloc (pp, pp);
		gsl_multifit_function_fdf ff;
		gsl_vector *x; 
		x = gsl_vector_alloc(pp);
		gsl_vector_set(x,0,f[i].tm3);
		gsl_vector_set(x,1,f[i].c);
		gsl_vector_set(x,2,-70);
		s = gsl_multifit_fdfsolver_alloc (T, tot_points, pp);
		// copy data to function
		ff.f = &eq_fit_trip;
		ff.df = NULL;
		ff.fdf = NULL;
		ff.p = pp;
		// total number of points is total of points
		// in the curve plus the number of points for the inv. fit
		ff.n = tot_points;
		ff.params = &d;
		gsl_multifit_fdfsolver_set (s, &ff, x);

		iter=0;
		do
		{
			iter++;
			status = gsl_multifit_fdfsolver_iterate (s);
			if(bVerbose){
				printf ("iter: %3u x = % 15.8f % 15.8f %15.8f "
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
		while (status == GSL_CONTINUE && iter < 500);
		gsl_multifit_covar (s->J, 0.0, covar);
		gsl_matrix_free (covar);
		gsl_vector_free(x);
	// copy tm3 data adjusted from each curve
		local_f[i].tm3 = gsl_vector_get(s->x, 0);
	}

	//free first solver
	gsl_multifit_fdfsolver_free (s);

	// do the 1/tm vs ln(ct) fitting
	const gsl_multifit_fdfsolver_type *Tl;
	gsl_multifit_fdfsolver *sl;
	// fit params in the straight line
	int ppl = 2; 
	gsl_matrix *covarl = gsl_matrix_alloc (ppl, ppl);
	struct fit_data dl = { n,i, local_f};
	gsl_multifit_function_fdf ffl;
	gsl_vector *xl;
	xl = gsl_vector_alloc(ppl);
	// DH and DS
	gsl_vector_set(xl,0,-70);
	gsl_vector_set(xl,1,-0.1);
	Tl = gsl_multifit_fdfsolver_lmsder;
	sl = gsl_multifit_fdfsolver_alloc (Tl, n, ppl);
	// copy data to function
	ffl.f=&eq_fit_straight;
	ffl.df = NULL;
	ffl.fdf = NULL;
	ffl.p = ppl;
	// total number of points the number of curves
	ffl.n = n;
	ffl.params = &dl;
	gsl_multifit_fdfsolver_set (sl, &ffl, xl);

	iter=0;
	do
	{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (sl);
		if(bVerbose){
			printf ("iter: %3u x = % 15.8f % 15.8f "
					"|f(x)| = %g\n",iter,
					gsl_vector_get (sl->x, 0),
					gsl_vector_get (sl->x, 1),
					gsl_blas_dnrm2 (sl->f));
		}

		if (status)
			break;
		status = gsl_multifit_test_delta (sl->dx, sl->x,
				1e-8, 1e-8);
	}
	while (status == GSL_CONTINUE && iter < 500);
	gsl_multifit_covar (sl->J, 0.0, covarl);

	#define FIT(i) gsl_vector_get(sl->x, i)
	#define ERR(i) sqrt(gsl_matrix_get(covarl,i,i))

	// compute contribution of inverse temperature to TSS
	for(i=0;i<n;i++){
		inv_t += ((real)1.0/(real)local_f[i].tm3);
		tot_inv_t++;
	}
	mean_inv_t = inv_t / (real)tot_inv_t;
	for(i=0;i<n;i++){
		tss_inv += (1.0/(real)local_f[i].tm3 - mean_inv_t ) * (1.0/(real)local_f[i].tm3 - mean_inv_t);
	}

	if (bUpdate){
		fit.dh3 = gsl_vector_get(sl->x, 0);
		fit.ds3 = gsl_vector_get(sl->x, 1);
		fit.dg3   = fit.dh3 - 298.15*fit.ds3;
		for(i=0; i<n; i++){
			f[i].tm3 = local_f[i].tm3 ;
		}
	}

	tss_tot = tss_inv ;

	double chi = gsl_blas_dnrm2(sl->f);
	fit.r2 = 1.0 - ( chi*chi / tss_tot ) ;
	double dof = n - ppl;
	double c = GSL_MAX_DBL(1, chi / sqrt(dof));

	if(bVerbose)
	{
		printf ("chisq/dof = %g\n",  pow(chi, 2.0) / dof);
		printf ("r2      = %g\n",  fit.r2);
		printf ("DH3    = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
		printf ("DS3    = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
		printf ("DG3    = %.5f +/- %.5f\n", FIT(0)-298*FIT(1),c*ERR(1)+c*298*ERR(0));
		printf ("status = %s\n", gsl_strerror (status));
	}

	gsl_multifit_fdfsolver_free (sl);
	gsl_matrix_free (covarl);
	gsl_vector_free(xl);

	return fit;
}

void print_multifit_info(f_info *f, double *m, double *h,  int n)
{
	int i;
	for(i=0; i<n; i++){
		printf("------------------ Curve %d ---------------------\n",i); 
		printf("Abs. max is %5.3f\n",m[i]); 
		printf("We suggest %5.3f as the triplex melting point\n",f[i].tm3); 
		printf("Transition height is %5.3f, ratio %5.3f \n",h[i], f[i].c); 
		printf("-------------------------------------------------\n"); 
	}
}
void scale_data(f_info *f,int n, double *m)
{
    int i;
    double min=0;
	for(i=0; i<n; i++){
		min = f[i].x[1][0];
		for(j=0; j<f[i].nts; j++){
			f[i].x[0][j] += 273.15;
			f[i].x[1][j] = (f[i].x[1][j] - min) / (m[i] - min)   ;
		}
	}

}

double *find_t_height(f_info *f,int n,double *m)
{
	int i,j; 
	double val=0;
	double *h;
	snew(h,n);
	for(i=0; i<n; i++){
		val=f[i].x[1][n-1];
		for(j=0; j<n; j++){
			if ( f[i].x[1][j] <= val ) { val = f[i].x[1][j] ; } 
		}
		h[i] = m[i] - val ;
	}
	return h;
}

void find_c_fix(f_info *f,int n,double *h)
{
	int i,j; 
	double val=0;
	double *cf;
	for(i=0; i<n; i++){
		f[i].c = (f[i].th2)/ (h[i] - f[i].th2 ) ;
		if (f[i].c<0) { 
			f[i].c=1; 
			printf("WARNING! Wrong with initial transition height\n");
			printf("WARNING! Setting the ratio of transitions to 1 to proceed.\n");
		} 
	}

}

double *find_maxabs(f_info *f, int n)
{
	int i,j;
	double val=0;
	double *v;
	snew(v,n);
	for(i=0; i<n; i++){
		for(j=0; j<f[i].nts; j++){
			if ( f[i].x[0][j] <= f[i].ef) { val = f[i].x[1][j] ; } 
		}
		v[i]=val; 	
	}
	return v; 
}

void compute_derivatives(double *x[], double *xp[], int nts)
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


void suggest_tm3_dh3(f_info *f,int nf, bool bVerbose)
{
	int i,j,k;
	double tm3=0;
	double val=0;
	double *dctrip[2];
	double *smth[2];
	// left/right number of frames to fullfill a local maxima
	int PAD=10;
	int lmax=0;
	double nts;
	int rmax=0;
	bool *bMax;

	for(k=0;k<nf;k++){
		// idk if they have been initialized to 0, but just in case
		// SHOULD I RESET ALL TO ZERO AT EACH ITERATION?
		nts=f[k].nts;
		snew(dctrip[0],nts-1);
		snew(dctrip[1],nts-1);
		snew(smth[0],nts-1);
		snew(smth[1],nts-1);
		snew(bMax,nts-1);
		for(i=0; i< nts-1; i++){ bMax[i]=FALSE; }
		compute_derivatives(f[k].x,dctrip,nts);
		// smooth a few times
		sgsmooth(dctrip,nts,smth); 
		sgsmooth(smth,nts,smth); 
		sgsmooth(smth,nts,smth); 
		sgsmooth(smth,nts,smth); 
		// search for local maxima within a window of 2*PAD
		// we addjust the padding based on the number of points
		PAD = (int) (nts / 10) ; 
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
		val=0;
		for(i=PAD; i< nts-PAD-1; i++){ 
			if (smth[0][i] < f[k].tm2-10 && bMax[i] && smth[1][i] > val) {
				val=smth[1][i];
				tm3=smth[0][i];
			}
		}
		if(tm3 == 0 || tm3 > f[k].tm2-5) {tm3=50+273.15;}
		f[k].tm3=tm3;

		sfree(dctrip[0]);
		sfree(dctrip[1]);
		sfree(smth[0]);
		sfree(smth[1]);
	}

}

f_all_sol opt_melt_diff_conc( f_info *f_inp, int nf,fit_option f_parm, bool bVerbose, const char *out )
{
	//
	// compute the ratio of transition heights and
	// the initial guess for tm3 
	// then fit all melting curves taking into account that
	// we want to have a straight line 1/Tm vs log(conc) 
	//
	int i=0;
	const char	*fitfile=NULL;
	const char	*out_file=NULL; 
	int file_counter=0;
	int iter=0;
	real *t_min;
	real *t_max;
	real *t_min_opt;
	real *t_max_opt;
	real r_number=0;
	real r_sq_min=0;
	real *ftriplex;
	FILE *f_fit;
	double *max_abs;
	double *t_height;
	f_all_sol mf;
	bool bUpdate;

	snew(ftriplex,nf);
	snew(t_min,nf);
	snew(t_max,nf);
	snew(t_min_opt,nf);
	snew(t_max_opt,nf);

	// remember to free max_abs and t_height and c_fix
	max_abs  = find_maxabs(f_inp,nf);
	t_height = find_t_height(f_inp,nf,max_abs);
	find_c_fix(f_inp,nf,t_height);
	//for(i=0;i<nf;i++){printf("%lf\n", c_fix[i]); } 
	scale_data(f_inp,nf,max_abs);
	suggest_tm3_dh3(f_inp,nf,bVerbose);
	if(bVerbose){
		print_multifit_info(f_inp, max_abs, t_height, nf);
	}
	sfree(t_height);
	sfree(max_abs);

	// Do not update the tm2 after each fitting attempt
	bUpdate = FALSE;
	for(iter=0; iter<f_parm.max_MC_iter; iter++){
		printf("\rProgress %3.2f", (real)iter/(real) f_parm.max_MC_iter);
		fflush(stdout);
		for(i=0; i<nf; i++){
			r_number =  (real)rand()/(real)RAND_MAX;
			t_min[i] = (f_inp[i].bf - f_parm.wt) + r_number*2*f_parm.wt ;
			r_number =  (real)rand()/(real)RAND_MAX;
			t_max[i] = (f_inp[i].ef - f_parm.wt) + r_number*2*f_parm.wt ;
		}
		mf=do_fit_triplex(f_inp, nf, t_min, t_max,  bVerbose, bUpdate);
		if ( mf.r2 > 0 && mf.r2 > r_sq_min ){
			for(i=0; i<nf; i++){
				t_min_opt[i] = t_min[i];
				t_max_opt[i] = t_max[i];
			}
			r_sq_min = mf.r2;
		}

	}

	// We should have found the optimals, keep the fit
	bUpdate = TRUE;
	mf=do_fit_triplex(f_inp, nf, t_min_opt, t_max_opt, bVerbose, bUpdate);
	eval_fit_straight(mf, f_inp, nf, ftriplex);
	// Output experimental vs fitted
	if(out) {
		f_fit = ffopen(out,"w");
		fprintf(f_fit,"# Opt fit DG %4.2f DH %4.2f DS %4.3f  -  R_2 %4.2f\n",mf.dg3, mf.dh3, mf.ds3, mf.r2) ;
		for(j=0; j<nf; j++){
			fprintf(f_fit,"%f %f \n", log(f_inp[j].conc), 1.0/f_inp[j].tm3);
		}

		fprintf(f_fit,"&\n");
		for(j=0; j<nf; j++){
			fprintf(f_fit,"%f %f \n", log(f_inp[j].conc), ftriplex[j]);
		}
		fclose(f_fit);
	}
	//Optim T solutions
	printf("Opt fit DG %4.2f DH %4.2f DS %4.3f  -  R_2 %4.2f\n",mf.dg3, mf.dh3, mf.ds3, mf.r2) ;
	if (bVerbose){
		printf("\n") ;
		for(i=0; i<nf; i++){
			printf(" %d - T_min %4.2f T_max %4.2f\n",i,t_min_opt[i], t_max_opt[i]) ;
		}
	}

	sfree(t_min);
	sfree(t_max);
	sfree(t_min_opt);
	sfree(t_max_opt);
	sfree(ftriplex);

	return mf;

}

int main(int argc, char *argv[])
{
	const char         *desc[] = {
		"Read triplex melting curve and duplex data to fit equilibrium [BR]"
			"Optimizes the best fit for 1/tm vs ln(Ct)."
	};

	FILE	*f_melt=NULL;
	FILE	*f_fit=NULL;
	static fit_option f_parm;
	static bool	bVerbose=FALSE;
	f_parm.max_MC_iter = 200;
	f_parm.wt = 2.0;

	/* Command-line arguments */
	t_pargs          pa[] = {
		{"-mi", FALSE, etINT, {&(f_parm.max_MC_iter)},
			"Number of iterations with random boundaries to find best r^2"},
		{"-wt", FALSE, etREAL, {&(f_parm.wt)},
			"Width of T window to randomize T_min/T_max boundaries"},
		{"-v", TRUE, etBOOL, {&bVerbose},
			"Be verbose."},
	};

	/* Output files */
	t_filenm            fnm[] = {
		{ efDAT, "-d", "triplex_input.dat", ffREAD },
		{ efXVG, "-o", "fitted_melt", ffWRITE  },
	};

#define NFILE asize(fnm)
	output_env_t    oenv;
	int				i,j,k;
	int				npargs;
	char			buf[256];
	const char		*inputfile=NULL;
	const char		*outfile=NULL;
	char			line[255];
	bool			bsettmax;
	f_info			*f_inp; 
	f_all_sol		all_sol;
	int				nfiles=0;
	int				nfiles_read=0;

	npargs=asize(pa);
	//	CopyRight(stderr, argv[0]);
	parse_common_args(&argc, argv, 0,
			NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);

	// read the input file
	inputfile = opt2fn("-d",NFILE,fnm);
	if(inputfile) {f_melt = ffopen(inputfile,"r");}
	while(fgets(line, 255, f_melt) != NULL) { nfiles++ ;}
	fclose(f_melt);
	snew(f_inp,nfiles);

	// reading all the inputs, f_inp contains the input file names, the 
	// number of lines per file, as well as the parameters for each fitting
	readfiles(inputfile, &f_inp, &nfiles_read);

	//initialize PRNG
	srand(time(NULL));

	// beautiful
	outfile = opt2fn("-o",NFILE,fnm);
	all_sol=opt_melt_diff_conc( f_inp, nfiles, f_parm, bVerbose, outfile );


	return 0;

}


