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
#include <sg.h>
#include <read_mult_t_inp.h>

#define R (0.0019858)
struct fit_params
{
    double k2, k3;
};

struct fit_data
{
    size_t n;
	f_info *f;
};

typedef struct{
	double dh3;
	double ds3;
	double dg3;
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
    int max_iter=100;
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
        //  if (status == GSL_SUCCESS)
        //  printf ("SS %d %f %f %f\n",i, a3, pf.k2,pf.k3);
        //  printf ("%5d [%.7f, %.7f] %.7f %.7f\n", iter, x_lo, x_hi, a2, x_hi-x_lo);
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

int eq_fit (const gsl_vector * x, void *d,
        gsl_vector *fit)
{
	int i;
	size_t n = ((struct fit_data *)d)->n;
	double dh3   = gsl_vector_get (x,2*n);
	double ds3   = gsl_vector_get (x,2*n+1);
	double *fd; // fraction duplex
	double *ft; // fraction triplex
	int ntot=0;
	for(i=0;i<n; i++){
		// fitting parameters
		double tm3 = gsl_vector_get (x,i);
		double c = gsl_vector_get (x,i+n);
		// constants and data to fit
		double tm2 = ((struct fit_data *)d)->f[i].tm2;
		double dh2 = ((struct fit_data *)d)->f[i].dh2;
		double nts = ((struct fit_data *)d)->f[i].nts;
		double *xx = ((struct fit_data *)d)->f[i].x[0];
		double *yy = ((struct fit_data *)d)->f[i].x[1];
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
			gsl_vector_set (fit, j+i, (Yi - yy[j]) );
			ntot++;
		}
		sfree(fd);
		sfree(ft);
	}
	for (i = 0; i < n; i++){
		double tm3 = gsl_vector_get (x,i);
		double conc = ((struct fit_data *)d)->f[i].conc;
		gsl_vector_set (fit,ntot+1, 
				(1.0/tm3  - ((R/dh3)*log(conc) + ((ds3-R*log(4))/dh3)) )  );
	}

    return GSL_SUCCESS;

}

f_all_sol do_fit_triplex(f_info *f,int n, bool bVerbose)
{
	int i;
	int iter=0;
	int tot_points=0;
	int status;
	f_all_sol fit;
	size_t pp = 0;
	// the number of parameters to fit is 2*n+DH+DS 
	// 2*n is for the tm3+c, 
	pp = 2*n+2;
	// Total number of points to fit
	for(i=0;i<n;i++){ tot_points+= f[i].nts; }
	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;

	gsl_matrix *covar = gsl_matrix_alloc (pp, pp);
	struct fit_data d = { n, f};
	gsl_multifit_function_fdf ff;

	gsl_vector *x; 
	x = gsl_vector_alloc(pp);
	for(i=0;i<n;i++){
		gsl_vector_set(x,i,f[i].tm3);
	}
	for(i=0;i<n;i++){
		gsl_vector_set(x,i+n,f[i].c);
	}
	// DH and DS
	gsl_vector_set(x,2*n,-70);
	gsl_vector_set(x,2*n+1,-0.1);

	T = gsl_multifit_fdfsolver_lmsder;
	s = gsl_multifit_fdfsolver_alloc (T, tot_points+n, pp);

	// copy data to function
	ff.f = &eq_fit;
	ff.df = NULL;
	ff.fdf = NULL;
	ff.p = pp;
	// total number of points is total of points
	// in the curve plus the number of points for the inv. fit
	ff.n = tot_points+n;
	ff.params = &d;

	gsl_multifit_fdfsolver_set (s, &ff, x);

	do
	{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);
		if(bVerbose){
			printf ("iter: %3u x = % 15.8f % 15.8f "
					"|f(x)| = %g\n",iter,
					gsl_vector_get (s->x, 2*n),
					gsl_vector_get (s->x, 2*n+1),
					gsl_blas_dnrm2 (s->f));
		}

		if (status)
			break;
		status = gsl_multifit_test_delta (s->dx, s->x,
				1e-4, 1e-4);
	}
	while (status == GSL_CONTINUE && iter < 500);
	gsl_multifit_covar (s->J, 0.0, covar);

#define FIT(i) gsl_vector_get(s->x, 2*n+i)
#define ERR(i) sqrt(gsl_matrix_get(covar,2*n+i,2*n+i))

	fit.dh3 = gsl_vector_get(s->x, 2*n);
	fit.ds3 = gsl_vector_get(s->x, 2*n+1);
	fit.dg3   = fit.dh3 - 298*fit.ds3;

	if(bVerbose)
	{
		double chi = gsl_blas_dnrm2(s->f);
		double dof = n - pp;
		double c = GSL_MAX_DBL(1, chi / sqrt(dof));

		printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);
		printf ("DH3    = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
		printf ("DS3    = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
		printf ("DG3    = %.5f +/- %.5f\n", fit.dg3,c*ERR(1)+c*298*ERR(0));
		printf ("status = %s\n", gsl_strerror (status));
	}

	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);

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

f_all_sol opt_melt_diff_conc( f_info *f_inp, int nf, bool bVerbose )
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
	double *max_abs;
	double *t_height;
	f_all_sol mf;

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

	//do the fit
	//	here you have to add the straight line fitting, besides working with all
	mf=do_fit_triplex(f_inp, nf, bVerbose);
	//	create_fit_curve(f_inp[i].x,ftrip,f_inp[i].nts,solution);

	/*
	// Output experimental vs fitted
	fitfile = opt2fn("-o",NFILE,fnm);
	char * s_fit = NULL;
	asprintf(&s_fit, "%s%d", fitfile, file_counter);
	if(s_fit) {
		f_fit = ffopen(s_fit,"w");
		for(j=0; j<f_inp[0].nts; j++){  
			fprintf(f_fit,"%f %f \n", f_inp[0].x[0][j], f_inp[0].x[1][j]);
		}
		fprintf(f_fit,"&\n");
		for(j=0; j<f_inp[0].nts; j++){  
			fprintf(f_fit,"%f %f \n", ftrip[0][j], ftrip[1][j]);
		}
		fclose(f_fit);
	}

	file_counter++;
	sfree(ftrip[0]);
	sfree(ftrip[1]);
	*/

	return mf;
}

int main(int argc, char *argv[])
{
	const char         *desc[] = {
		"Read triplex melting curve and duplex data to fit equilibrium [BR]"
			"Solving the equations takes time."
	};

	FILE	*f_melt;
	FILE	*f_fit;
	static bool	bVerbose=FALSE;

	/* Command-line arguments */
	t_pargs          pa[] = {
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
	readfiles(inputfile, f_inp, &nfiles_read);

	if( nfiles_read != nfiles ) {gmx_fatal(FARGS,"Wrong number of lines read\n"); }

	// beautiful
	all_sol=opt_melt_diff_conc( f_inp, nfiles, bVerbose );


	return 0;

}


