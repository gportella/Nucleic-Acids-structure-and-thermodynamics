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
#include "read_mult_d_inp.h"

#define R (0.0019858)
struct fit_params
{
    double k2, k3;
};

struct fit_data
{
    size_t n;
    int  w;
	f_info *f;
};

typedef struct{
	double dh2;
	double ds2;
	double dg2;
	double r2;
} f_all_sol; 

typedef struct{
	int max_MC_iter;
	real wt;
} fit_params;

typedef struct
{
	double tm2;
	double dh2;
	double tm3;
	double dh3;
	double c;
	double dg3;
} f_p;

void do_solve_eq_d(double *x, int n, f_p p, double *fd);


void create_fit_curve(double *x[],double *xf[],int n,f_p s)
{
    int i;
    double *fd; // fraction duplex
    snew(fd,n);
    do_solve_eq_d(x[0],n,s,fd);
    for (i = 0; i < n; i++)
    {
        //  ( 1 - a2 + c * ( 1- a3) )/(1+c)
        xf[1][i] = ( 1 - fd[i] + s.c * ( 1- fd[i]) )/(1+s.c);
        xf[0][i] = x[0][i];
    }

}

void do_solve_eq_d(double *x, int n, f_p p, double *fd)
{
	double *k;
	snew(k,n);
	int max_iter=100;
	for(i=0;i<n; i++){
		k[i] = exp((p.dh2/R)*(1.0/p.tm2 - 1.0/x[i]));
		fd[i] = (sqrt( 1.0 + 2*k[i])-1)/(k[i]);
	}
	sfree(k);
}

int eq_fit (const gsl_vector * x, void *d,
        gsl_vector *fit)
{
	int i;
	size_t n = ((struct fit_data *)d)->n;
	int  w = ((struct fit_data *)d)->w;
	double tm2 = gsl_vector_get (x,0);
	double c = gsl_vector_get (x,1);
	double dh2   = gsl_vector_get (x,2);
	double ds2   = gsl_vector_get (x,3);
	double *fd; // fraction duplex
	// fitting parameters
	// constants and data to fit
	double nts = ((struct fit_data *)d)->f[w].nts - 1 ;
	double *xx = ((struct fit_data *)d)->f[w].xp[0];
	double *yy = ((struct fit_data *)d)->f[w].xp[1];
	f_p fpt;
	fpt.tm2=tm2;
	fpt.dh2=dh2;
	snew(fd,nts);

	do_solve_eq_d(xx,nts,fpt,fd);

	for (j = 0; j < nts; j++)
	{
		//  acc. to R&C PNAS 1991 it should have a diff 
		//  form, but this one from a triplex paper fits better
		//  ( 1 - a2 + c * ( 1- a3) )/(1+c)
		if ( c < 0 ) {c=100000;}
		double Yi = (  c*(((1-fd[j])*fd[j])/(2-fd[j]))*xx[j]*xx[j] );
		gsl_vector_set (fit, j+i, (Yi - yy[j]) );
	}

	sfree(fd);

    return GSL_SUCCESS;

}

// function that the gls solver calls to make the fit
// for 1/tm2 vs ln(Ct)l 
int eq_fit_straight(const gsl_vector * x, void *d,
        gsl_vector *fit)
{
    int i;
    double dh2 = gsl_vector_get (x,0);
    double ds2 = gsl_vector_get (x,1);
    size_t n = ((struct fit_data *)d)->n;
    for (i = 0; i < n; i++){
        double conc = ((struct fit_data *)d)->f[i].conc;
        double tm2 = ((struct fit_data *)d)->f[i].tm2;
        gsl_vector_set (fit,i,
                (1.0/tm2  - ((R/dh2)*log(conc) + ((ds2-R*log(4))/dh2)) )  );
    }

    return GSL_SUCCESS;
}


void eval_fit_straight(f_all_sol mf, f_info *f_inp, int nf, real *f)
{
	int i; 
	for(j=0; j<nf; j++){
		f[j] = (R/mf.dh2)*log(f_inp[j].conc) + ((mf.ds2-R*log(4))/mf.dh2);
	}

}

f_all_sol do_fit_duplex(f_info *f, int n, real *t_min, real *t_max, bool bVerbose, bool bUpdate)
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

	snew(local_f,n);

	// count points within boundaries, allocate 
	// and copy to new array. Only copy the xp (derivative)
	// as is the one that matters for the fitting
	for(i=0; i<n; i++){
		for(j=0; j<f[i].nts; j++){
			if( f[i].xp[0][j] >= t_min[i] && f[i].xp[0][j] <= t_max[i] ) { local_f[i].nts++ ; }
		}
		snew(local_f[i].xp[0],local_f[i].nts); 
		snew(local_f[i].xp[1],local_f[i].nts); 
	}
	// now copy the data
	// not only the derivatives of the absorbance, also the conc
	k=0;
	for(i=0; i<n; i++){
		local_f[i].conc = f[i].conc ;
		for(j=0; j<f[i].nts; j++){
			if( f[i].xp[0][j] >= t_min[i] && f[i].xp[0][j] <= t_max[i] ) { 
				local_f[i].xp[0][k] = f[i].xp[0][j] ;
				local_f[i].xp[1][k] = f[i].xp[1][j] ;
				// we use this loop to compute the TSS, the total sum of squares of the
				// "y" data, to be used for r_squared after we know chi_sq
				abs_t += local_f[i].xp[1][k] ;
				tot_abs++;
				k++;
			}
		}
		k=0;
	}

	// Total number of points to fit
	size_t pp = 4;
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
        gsl_vector_set(x,0,f[i].tm2);
        gsl_vector_set(x,1,f[i].c);
        gsl_vector_set(x,2,-70);
        gsl_vector_set(x,3,-0.1);
        s = gsl_multifit_fdfsolver_alloc (T, tot_points, pp);
        // copy data to function
        ff.f = &eq_fit;
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
    // copy tm2 data adjusted from each curve
        local_f[i].tm2 = gsl_vector_get(s->x, 0);
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
        inv_t += ((real)1.0/(real)local_f[i].tm2);
        tot_inv_t++;
    }
    mean_inv_t = inv_t / (real)tot_inv_t;
    for(i=0;i<n;i++){
        tss_inv += (1.0/(real)local_f[i].tm2 - mean_inv_t ) * (1.0/(real)local_f[i].tm2 - mean_inv_t);
    }

    if (bUpdate){
        fit.dh2 = gsl_vector_get(sl->x, 0);
        fit.ds2 = gsl_vector_get(sl->x, 1);
        fit.dg2   = fit.dh2 - 298.15*fit.ds2;
        for(i=0; i<n; i++){
            f[i].tm2 = local_f[i].tm2 ;
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

void print_multifit_info(f_info *f, int n)
{
	int i;
	for(i=0; i<n; i++){
		printf("------------------ Curve %d ---------------------\n",i); 
		printf("We suggest %5.3f as the duplex melting point\n",f[i].tm2); 
		printf("-------------------------------------------------\n"); 
	}
}
void scale_data(f_info *f,int n)
{
    int i,j;
	for(i=0; i<n; i++){
		for(j=0; j<f[i].nts; j++){
			f[i].x[0][j] += 273.15;
		}
	}

}


void find_c_fix(f_info *f,int n)
{
	int i,j; 
	for(i=0; i<n; i++){
		f[i].c = 0.0001 ;
	}

}


void compute_derivatives(f_info *f, int nc)
{
	int i=0;
	int j=0;
	double x0=0;
	double y0=0;

	for(j=0; j<nc; j++){
		x0=f[j].x[0][0];
		y0=f[j].x[1][0];
			for (i=1; i<f[j].nts; i++){
			f[j].xp[0][i-1] = x0+(f[j].x[0][i] -x0)/2 ; 
			f[j].xp[1][i-1] = (f[j].x[1][i] -y0)/(f[j].x[0][i]-x0) ;
			x0=f[j].x[0][i] ; 
			y0=f[j].x[1][i] ;
		}
	}
}


void suggest_tm2_dh2(f_info *f,int nf, bool bVerbose)
{
	int i,j,k;
	double tm2=0;
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
		snew(smth[0],nts-1);
		snew(smth[1],nts-1);
		snew(bMax,nts-1);
		for(i=0; i< nts-1; i++){ bMax[i]=FALSE; }
		// smooth a few times
		// I should probably smooth the data not the derivative
		// see what i did in python
		sgsmooth(f[k].xp,nts-1,smth); 
		sgsmooth(smth,nts-1,smth); 
		//sgsmooth(smth,nts-1,smth); 
		//sgsmooth(smth,nts-1,smth); 
		//for(i=0; i< nts-1; i++){ printf("%lf %lf \n", smth[0][i], smth[1][i]); }
		//for(i=0; i< nts-1; i++){ printf("%lf %lf \n", f[k].xp[0][i], f[k].xp[1][i]); }
		//printf("\n");
		// search for local maxima within a window of 2*PAD
		// we addjust the padding based on the number of points
		PAD = (int) (nts / 20) ; 
		//#pragma omp parallel for private(j,lmax,rmax)
		for(i=PAD; i< nts-PAD-1; i++){
			for(j=1; j<PAD; j++){
				if (smth[1][i-j] < smth[1][i] ) { lmax++; }
				if (smth[1][i+j] < smth[1][i] ) { rmax++; } 
			}
			if(lmax+rmax==2*(PAD-1)) { bMax[i] = TRUE ; }
			lmax=0;
			rmax=0;
		}
		//for(i=0; i< nts-1; i++){ printf("%f %f \n",smth[0][i],(float)0.001*bMax[i]);}
		// throw away after the duplex melting curve, keep whatever is first
		tm2=0;
		val=0;
		for(i=PAD; i< nts-PAD-1; i++){ 
			if (smth[0][i] < f[k].ef && bMax[i] && smth[1][i] > val) {
				val=smth[1][i];
				tm2=smth[0][i];
			}
		}
		if(tm2 == 0) {tm2=50+273.15;}
		f[k].tm2=tm2;

		sfree(smth[0]);
		sfree(smth[1]);
	}

}

f_all_sol opt_melt_diff_conc(f_info *f_inp, int nf, fit_params f_parm, bool bVerbose ,const char *out )
{
	//
	// Scale data 
	// fit all melting curves taking into account that
	// we want to have a straight line 1/Tm vs log(conc) 
	//
	int i=0;
	int iter=0;
	real *t_min;
	real *t_max;
	real *t_min_opt;
	real *t_max_opt;
	real r_number=0;
	real r_sq_min=0;
	FILE *f_fit;
	f_all_sol mf;
	real *fduplex;
	bool bUpdate;

	snew(fduplex,nf);
	snew(t_min,nf);
	snew(t_max,nf);
	snew(t_min_opt,nf);
	snew(t_max_opt,nf);

	// remember to free max_abs and t_height and c_fix
	find_c_fix(f_inp,nf);
	scale_data(f_inp,nf);
	compute_derivatives(f_inp,nf);
	suggest_tm2_dh2(f_inp,nf,bVerbose);
	if(bVerbose){
		print_multifit_info(f_inp, nf);
	}

	// we do a number of iterations
	// and find the boundary conditions
	// that give the best r_2

	// We do not update the tm2 guess in order to have
	// reproduceble fitting results. Once we have found
	// the best fit, we can keep the results
	bUpdate = FALSE;
	for(iter=0; iter<f_parm.max_MC_iter; iter++){
	//	printf("\rProgress %3.2f", (real)iter/(real) f_parm.max_MC_iter);
	//	fflush(stdout);
		for(i=0; i<nf; i++){
			r_number =  (real)rand()/(real)RAND_MAX;
			t_min[i] = (f_inp[i].bf - f_parm.wt) + r_number*2*f_parm.wt ;
			r_number =  (real)rand()/(real)RAND_MAX;
			t_max[i] = (f_inp[i].ef - f_parm.wt) + r_number*2*f_parm.wt ;
		}
		// Do the fit
		// Here you have to add the straight line fitting, besides working with all
		// It is here that we select based on the tmin/tmax 
		mf=do_fit_duplex(f_inp, nf, t_min, t_max, bVerbose, bUpdate);
		if ( mf.r2 > 0 && mf.r2 > r_sq_min ){
			for(i=0; i<nf; i++){
				t_min_opt[i] = t_min[i];
				t_max_opt[i] = t_max[i];
			}
			r_sq_min = mf.r2;
		}
	}
	printf("\n");

	// We should have found the optimals, keep the fit
	bUpdate = TRUE;
	mf=do_fit_duplex(f_inp, nf, t_min_opt, t_max_opt, bVerbose, bUpdate);
	eval_fit_straight(mf, f_inp, nf, fduplex);
	//create_fit_curve(f_inp[i].x,ftrip,f_inp[i].nts,solution);
	// Output experimental vs fitted
	if(out) {
		f_fit = ffopen(out,"w");
		for(j=0; j<nf; j++){  
			fprintf(f_fit,"%f %f \n", log(f_inp[j].conc), 1.0/f_inp[j].tm2);
		}

		fprintf(f_fit,"&\n");
		for(j=0; j<nf; j++){  
			fprintf(f_fit,"%f %f \n", log(f_inp[j].conc), fduplex[j]);
		}
		fclose(f_fit);
	}
	//Optim T solutions
	printf("Opt fit DG %4.2f DH %4.2f DS %4.2f  -  R_2 %4.2f\n",mf.dg2, mf.dh2, mf.ds2, mf.r2) ;
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
	return mf;
}

int main(int argc, char *argv[])
{
	const char         *desc[] = {
		"Read duplex melting curve and duplex data to fit equilibrium [BR]"
			"Solving the equations takes time."
	};

	FILE	*f_melt;
	FILE	*f_fit;
	static real	tmaxabs=-1;
	static bool	bVerbose=FALSE;
	static fit_params f_parm;
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
		{ efDAT, "-d", "duplex_input.dat", ffREAD },
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
	real 			*t_min;
	real 			*t_max;

	npargs=asize(pa);
	//	CopyRight(stderr, argv[0]);
	parse_common_args(&argc, argv, 0,
			NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);

	// read the input file
	inputfile = opt2fn("-d",NFILE,fnm);

	// reading all the inputs, f_inp contains the input file names, the 
	// number of lines per file, as well as the parameters for each fitting
	// We are going to read all the points here, then apply the
	// limits later on 
	readfiles(inputfile, &f_inp, &nfiles_read);

	//initialize PRNG
	srand(time(NULL));

	outfile = opt2fn("-o",NFILE,fnm);
	all_sol=opt_melt_diff_conc(f_inp, nfiles_read, f_parm,  bVerbose, outfile);

	return 0;

}


