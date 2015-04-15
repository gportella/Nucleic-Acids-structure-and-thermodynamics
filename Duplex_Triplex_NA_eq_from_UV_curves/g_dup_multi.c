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
#include <read_mult_d_inp.h>

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
	double dh2;
	double ds2;
	double dg2;
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
	double dh2   = gsl_vector_get (x,2*n);
	double ds2   = gsl_vector_get (x,2*n+1);
	double *fd; // fraction duplex
	int ntot=0;
	for(i=0;i<n; i++){
		// fitting parameters
		double tm2 = gsl_vector_get (x,i);
		double c = gsl_vector_get (x,i+n);
		// constants and data to fit
		double nts = ((struct fit_data *)d)->f[i].nts - 1 ;
		double *xx = ((struct fit_data *)d)->f[i].xp[0];
		double *yy = ((struct fit_data *)d)->f[i].xp[1];
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
			ntot++;
		}
		sfree(fd);
	}
	for (i = 0; i < n; i++){
		double tm2 = gsl_vector_get (x,i);
		double conc = ((struct fit_data *)d)->f[i].conc;
		gsl_vector_set (fit,ntot+1, 
				(1.0/tm2  - ((R/dh2)*log(conc) + ((ds2-R*log(4))/dh2)) )  );
	}

    return GSL_SUCCESS;

}

f_all_sol do_fit_duplex(f_info *f,int n, bool bVerbose)
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
		gsl_vector_set(x,i,f[i].tm2);
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

	fit.dh2 = gsl_vector_get(s->x, 2*n);
	fit.ds2 = gsl_vector_get(s->x, 2*n+1);
	fit.dg2   = fit.dh2 - 298*fit.ds2;

	for(i=0;i<n;i++){
		f[i].tm2 = gsl_vector_get(s->x, i);
	}

	if(bVerbose)
	{
		double chi = gsl_blas_dnrm2(s->f);
		double dof = n - pp;
		double c = GSL_MAX_DBL(1, chi / sqrt(dof));

		printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);
		printf ("DH3    = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
		printf ("DS3    = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
		printf ("DG3    = %.5f +/- %.5f\n", fit.dg2,c*ERR(1)+c*298*ERR(0));
		printf ("status = %s\n", gsl_strerror (status));
	}

	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);

	return fit;
}

void print_multifit_info(f_info *f, int n)
{
	int i;
	for(i=0; i<n; i++){
		printf("------------------ Curve %d ---------------------\n",i); 
		printf("We suggest %5.3f as the triplex melting point\n",f[i].tm2); 
		printf("-------------------------------------------------\n"); 
	}
}
void scale_data(f_info *f,int n)
{
    int i;
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

f_all_sol opt_melt_diff_conc( f_info *f_inp, int nf, bool bVerbose , const char *out )
{
//
// compute the ratio of transition heights and
// the initial guess for tm3 
// then fit all melting curves taking into account that
// we want to have a straight line 1/Tm vs log(conc) 
//
	int i=0;
	int file_counter=0;
	FILE *f_fit;
	f_all_sol mf;

	// remember to free max_abs and t_height and c_fix
	find_c_fix(f_inp,nf);
	//for(i=0;i<nf;i++){printf("%lf\n", c_fix[i]); } 
	scale_data(f_inp,nf);
	compute_derivatives(f_inp,nf);
	suggest_tm2_dh2(f_inp,nf,bVerbose);
	if(bVerbose){
		print_multifit_info(f_inp, nf);
	}

	//do the fit
	//	here you have to add the straight line fitting, besides working with all
	mf=do_fit_duplex(f_inp, nf, bVerbose);
	//	create_fit_curve(f_inp[i].x,ftrip,f_inp[i].nts,solution);

	// Output experimental vs fitted
	if(out) {
		f_fit = ffopen(out,"w");
		for(j=0; j<nf; j++){  
			fprintf(f_fit,"%f %f \n", log(f_inp[j].conc), 1.0/f_inp[j].tm2);
		}
		/*
		fprintf(f_fit,"&\n");
		for(j=0; j<f_inp[0].nts; j++){  
			fprintf(f_fit,"%f %f \n", ftrip[0][j], ftrip[1][j]);
		}
		*/
		fclose(f_fit);
	}

	file_counter++;

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

	/* Command-line arguments */
	t_pargs          pa[] = {
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

	outfile = opt2fn("-o",NFILE,fnm);
	all_sol=opt_melt_diff_conc( f_inp, nfiles, bVerbose, outfile );


	return 0;

}


