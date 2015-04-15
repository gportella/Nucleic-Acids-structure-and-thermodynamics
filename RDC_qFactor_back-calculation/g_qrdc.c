// GPC 10/03/2015 g_qrdc v 0.2 
//
// Code to back-compute RDCs the theta-method way. Finds the optimal correlation
// for computed vs exp RDCs by x/y rotations, then back-calculates the RDCs in
// the best orientation. Based on RDC colvar from plumed 
//
//
// Easily fixable caveats : No pbc taken into account, whole molecules only
// "Feature" : Only one scaling for all type of bonds
//
// TODO: 
// -- Add option not to redo aligment, needs to ask for scaling factor
// -- Better parsing and checking of the input RDC file
// -- Add more error checks for routines that could fail
// -- Add average + stdev for correlations found
// -- Add switch to be able to output or not the aligned xtc 
//
// RDC input file free format, but no comments or extra lines allowed. The
// format is 
//
// atom_index atom_index RDC gyormagnetic_ratio
// 
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

#define Const (0.3356806)

typedef struct{
	int ind[2];		// index for each atom in the bond / atom pair
	real exp;		// experimental RDC
	real *comp;		// computed RDC
	real r_comp;	// replica-averaged RDC
	real av_comp;	// trajectory-averaged RDC
	real gyr;		// Gyromagnetic ratio for the bond / atom pair
	real scale;		// Scaling factor between computed and experimental RDC
} rdc_str_e;


// copied from http://tinyurl.com/nh9xtas
// removed non-C99 assert call and placed an if instead.
int randint(int n) {
  if ((n - 1) == RAND_MAX) {
    return rand();
  } else {
    // Chop off all of the values that would cause skew...
    long end = RAND_MAX / n; // truncate skew
    if ((end > 0L)==FALSE) {printf("Sth terribly wrong happens with RAND MAX\n"); exit(1);}
    end *= n;

    // ... and ignore results from rand() that fall above that limit.
    // (Worst case the loop condition should succeed 50% of the time,
    // so we can expect to bail out of this loop pretty quickly.)
    int r;
    while ((r = rand()) >= end);
    return r % n;
  }
}

static void rotate(int natom, rvec *x, rvec *new_x, real alfa, real beta, real gamma)
{
//
// Rotate the rvec x about x/y/z angles. For our case we could remove the
// gamma, but anyway. Shamelessly copied from gmx. Added copy of original
// coordinates
//
    int  i;
    rvec x_old;

    for (i = 0; i < natom; i++)
    {
        copy_rvec(x[i], x_old);
        /*calculate new x[i] by rotation alfa around the x-axis*/
        new_x[i][XX] =   x_old[XX];
        new_x[i][YY] =             cos(alfa)*x_old[YY] - sin(alfa)*x_old[ZZ];
        new_x[i][ZZ] =             sin(alfa)*x_old[YY] + cos(alfa)*x_old[ZZ];
        copy_rvec(new_x[i], x_old);
        /*calculate new x[i] by rotation beta around the y-axis*/
        new_x[i][XX] =   cos(beta)*x_old[XX]           + sin(beta)*x_old[ZZ];
        new_x[i][YY] =                       x_old[YY];
        new_x[i][ZZ] = -sin(beta)*x_old[XX]           + cos(beta)*x_old[ZZ];
        copy_rvec(new_x[i], x_old);
        /*calculate new x[i] by rotation gamma around the z-axis*/
        new_x[i][XX] = x_old[XX]*cos(gamma) - x_old[YY]*sin(gamma);
        new_x[i][YY] = x_old[XX]*sin(gamma) + x_old[YY]*cos(gamma);
        new_x[i][ZZ] =                                             x_old[ZZ];
    }
}

real compute_av_qfact( rdc_str_e *rdc, int n_rdc, int nframes )
{
//
// Computes normalization of RDC averages by the number of frames and returns
// the Q-factor of the averged computed RDC 
//

	int i=0;
	real sum=0;
	real Q=0;
	real q_fact=0;

	for (i=0; i< n_rdc; i++) {
		rdc[i].av_comp /= nframes;
		sum+=(rdc[i].exp-rdc[i].av_comp)*(rdc[i].exp-rdc[i].av_comp);
		Q+=(rdc[i].exp*rdc[i].exp);
	}
	q_fact = sqrt(sum/Q); 
	return q_fact;
}
real compute_rdc_qfact(rdc_str_e *rdc, int n_rep, int n_rdc, int natoms, rvec **x, bool bCorrel)
{
//
// Returns either the correlation between computed and experimental RDC or the
// Q-factor between them. Recomended use is with scaling factors of 1 for
// bCorrel == True, and with the actual scaling factor otherwise.  
// Addapted from plumed's colvar RDC.
//
    int a1=0;
    int a2=0;
    int i=0;
    int r=0;
    real scx=0;
    real scx2=0;
    real scy=0;
    real scy2=0;
    real scxy=0;
    real sum=0;
    real Q=0;
    real corr_qfact=0;
    real r_norm=0;

    r_norm = 1. / n_rep;
    //we set the ensemble averaged RDC to zero 
    for (i=0; i< n_rdc; i++) { rdc[i].r_comp = 0; }

    for (r=0; r< n_rep; r++){
        for (i=0; i< n_rdc ; i++){
            // gmx counts from 0, the input should be compatible
            // with the given pdb file
            a1 = rdc[i].ind[0]-1 ;
            a2 = rdc[i].ind[1]-1 ;
            // this could be changed to account for PBCs 
            real xx =  x[r][a1][XX]-x[r][a2][XX] ;
            real yy =  x[r][a1][YY]-x[r][a2][YY] ;
            real zz =  x[r][a1][ZZ]-x[r][a2][ZZ] ;
            real d = sqrt( xx*xx + yy*yy + zz*zz );
            real d2   = d*d;
            real d3   = d2*d;
            real id3  = 1./d3;
            real max  = -Const*rdc[i].scale*rdc[i].gyr;
            real dmax = id3*max;
            real cos_theta = zz / d ;
            rdc[i].comp[r] = 0.5*dmax*(3.*cos_theta*cos_theta-1.);
			//printf("Exp is %f cal %f scal %f\n", rdc[i].exp, rdc[i].comp[r], rdc[i].scale);
            // replica averaged RDC
            rdc[i].r_comp += rdc[i].comp[r]*r_norm ;
        }
    }
    if (bCorrel){
        for (i=0; i< n_rdc ; i++){
            scx  += rdc[i].r_comp;
            scx2  += rdc[i].r_comp*rdc[i].r_comp;
            scy  += rdc[i].exp;
            scy2  += rdc[i].exp*rdc[i].exp;
            scxy += rdc[i].exp*rdc[i].r_comp;
        }
        real num = (n_rdc*n_rep)*scxy - scx*scy;
        real idevx = 1./sqrt((n_rdc*n_rep)*scx2-scx*scx);
        real idevy = 1./sqrt((n_rdc*n_rep)*scy2-scy*scy);
        corr_qfact = num * idevx * idevy;
    } else {
        for (i=0; i< n_rdc ; i++){
            // Add up the computed RDCs to get ensemble value later
            // WARNING!! make sure to CALL this ONCE PER FRAME 
            // or the STATISTICS WILL BE BOTCHED!!. 
            sum+=(rdc[i].exp-rdc[i].r_comp)*(rdc[i].exp-rdc[i].r_comp);
            Q+=(rdc[i].exp*rdc[i].exp);
            rdc[i].av_comp += rdc[i].r_comp;
        }
        // per frame ensemble averaged Q-factor
        corr_qfact = sqrt(sum/Q);
    }
    return corr_qfact;
}

void find_best_orient(rdc_str_e *rdc, int n_rep, int n_rdc, int natoms, rvec *x, int max_iter, 
		real mc_half_step, real inv_T, bool bVerbose)
{
//
// Find the alignment that best correlates the computed RDCs and the
// experimental ones. It uses a Metropolis Monte Carlo search. Parameters are
// tunable, but defaults should be okeish. After the alignment is found,
// computes the slope of the correlation (aka scaling) and computes the
// Q-factor. The rdc_str_e structure containts the best rdcs per frame, and
// a cumulative sum, such that we can perform the average later on. 
//
	int xrot=0;
	int count=0;
	int i=0;
	int j=0;
	int r=0;
	int r_rot=0;
	real corr=0;
	real corr_b4=0;
	real max_corr=0;
	real rotx=0;
	real roty=0;
	real xr=0;
	real yr=0;
	real ran_met=0;
	real slope=0;
	real *fx;
	real *fy;
	real qfact=0;
	rvec **xens;
	rvec *xnew;
	rvec **xnew_e;
	rvec **xopt_e;
	bool bCorrel=TRUE;

	//initiate PNRG
	srand(time(NULL));

	snew(xnew,natoms);
	snew(xopt_e,n_rep);
	snew(xens,n_rep);
	snew(xnew_e,n_rep);

	// create initial replicas of the frame
	for(r=0; r<n_rep; r++){
		snew(xens[r],natoms);	
		snew(xnew_e[r],natoms);	
		snew(xopt_e[r],natoms);	
		for(j=0; j<natoms; j++){
			copy_rvec(x[j],xens[r][j]);
			copy_rvec(x[j],xnew_e[r][j]);
			copy_rvec(x[j],xopt_e[r][j]);
		}
	}
	
	// set the scaling to 1 for correlation
	for (i=0; i< n_rdc; i++) { rdc[i].scale = 1  ; }
	// initiate qfact
	corr_b4 = compute_rdc_qfact(rdc, n_rep, n_rdc, natoms, xens, bCorrel);
	max_corr = corr_b4;
	int acc=0;
	do {
		count++;
		//select which one you want to rotate 
		r_rot = randint(n_rep) ;
		// generate random rotation [M,N] range
		// M + rand() / (RAND_MAX / (N - M + 1) + 1)
		xr = -mc_half_step + ((real) rand() / ( (real) RAND_MAX / (mc_half_step*2)+1 )) ;
		yr = -mc_half_step + ((real) rand() / ( (real) RAND_MAX / (mc_half_step*2)+1 )) ;
		rotx = DEG2RAD * xr;
		roty = DEG2RAD * yr;
		// propose new rotation
		rotate(natoms, xens[r_rot], xnew, rotx, roty, 0);
		//copy xnew to its place in xnew_e
		for(i=0; i<natoms; i++){ copy_rvec(xnew[i],xnew_e[r_rot][i]); }
		// Try qfact rotated ala Metropolis
		corr = compute_rdc_qfact(rdc, n_rep, n_rdc, natoms, xnew_e, bCorrel);
		ran_met = ((real) rand() / (real) RAND_MAX )*(1);
		if (exp((corr_b4-corr)/inv_T) < ran_met  ){
			if (corr > max_corr){
				max_corr = corr ; 
				// We keep it if it improves
				// I don't do anything with them at the moment
				for(i=0;i<natoms; i++){copy_rvec(xnew_e[r_rot][i],xopt_e[r_rot][i]);}
			}
			// Copy it for next iteration
			for(i=0;i<natoms; i++){copy_rvec(xnew_e[r_rot][i],xens[r_rot][i]);}
			corr_b4 = corr ; 
			acc++;
		}
	} while ( max_corr < 0.999  && count < max_iter );

	// rdc structure has been updated with newest computed values
	corr=compute_rdc_qfact(rdc, n_rep, n_rdc, natoms, xopt_e, bCorrel);
	if (bVerbose) {
		printf ("\nMax. correlation after MMC is %f (%d iterations)\n", max_corr, count);
	} 
	// fill up arrays for fitting
    snew(fx,n_rdc);
    snew(fy,n_rdc);
	for(i=0; i<n_rdc;i++){
		fx[i]=rdc[i].r_comp;
		fy[i]=rdc[i].exp;
	}
    lsq_y_ax( n_rdc, fx, fy, &slope); // fitting  y=ax
	sfree(fx);
	sfree(fy);
	// the slope is the scaling factor
	for (i=0; i< n_rdc; i++) { rdc[i].scale = slope ; }
	// Now we do want to compute a Q-factor
	// The rdc addition for average calculation takes place
	bCorrel=FALSE;
	qfact=compute_rdc_qfact(rdc, n_rep, n_rdc, natoms, xopt_e, bCorrel);
	if (bVerbose){
		printf("Per frame Q-fact is %f \n",qfact);
	}
	// Here we could copy new orientation for output gmx traj
	// This would require n_rep output traj.
	sfree(xnew);
	for(r=0; r< n_rep; r++){
		sfree(xnew_e[r]);
		sfree(xopt_e[r]);
		sfree(xens[r]);	
	}

}

void write_backrdc(FILE *out, rdc_str_e *rdc,int n_rdc, real av_qfact)
{
	int i=0;
	
	fprintf(out, "# Q-factor %6.3f\n", av_qfact);
	for (i=0; i< n_rdc; i++){
		fprintf(out, "%f  %f\n", rdc[i].av_comp, rdc[i].exp);
	}

}


int main(int argc, char *argv[])
{
	const char         *desc[] = {
		"Read a pdb + trajectory, a set of RDCs and compute the Q-factor[BR]"
			"between the ensemble of structures and the experimental values."
	};

	FILE				*f_rdc;
	FILE				*f_back;
	static int	max_iter=1000;
	static int	n_rep=4;
	static real	mc_half_step=5;
	static real	inv_T=0.0001;
	static bool	bVerbose=FALSE;

	/* Command-line arguments */
	t_pargs          pa[] = {
		//TODO: implement that later on
	//	{"-nomc", FALSE, etBOOL, {&bMC},
	//		"Do not redo alignment"},
		{"-mi", FALSE, etINT, {&max_iter},
			"Maximum number of iterations MC search"},
		{"-ms", FALSE, etREAL, {&mc_half_step},
			"MC half step size, in degrees"},
		{"-mt", FALSE, etREAL, {&inv_T},
			"MC weight factor / inverse temperature"},
		{"-nr", FALSE, etINT, {&n_rep},
			"Number of replicas"},
		{"-v", TRUE, etBOOL, {&bVerbose},
			"Be verbose, output maximum correlation for each frame to std output."},
	};

	/* Output files */
	t_filenm            fnm[] = {
		{ efTRX, "-f", NULL, ffREAD },
		{ efTPS, NULL, NULL, ffREAD },
		{ efDAT, "-r", "rdc.dat", ffREAD },
		//{ efTRX, "-o", "trajout", ffWRITE  },
		{ efXVG, "-q", "comp_rdc", ffWRITE  },
	};

#define NFILE asize(fnm)
	output_env_t	oenv;
	t_tpxheader		header;
	matrix			box;
	t_topology      top;
	t_trxstatus		*status;
	t_fileio		*trxout;
	real            t;
	rvec			*x,*xp;
	real			av_qfact=0;
	int				natoms,i,j,kkk;
	int				ePBC;
	int				nframes=0;
	int				n_rdc=0;
	int				step=1;
	char			buf[256];
	const char		*rdcfile=NULL;
	const char		*backfile=NULL;
	const char		*out_file=NULL; 
	char			line[255];


	// The structure that will contain the RDC info
	// check definition at the begining
	rdc_str_e *rdc;


	int npargs;
	npargs=asize(pa);
	CopyRight(stderr, argv[0]);
	parse_common_args(&argc, argv, PCA_CAN_TIME | PCA_BE_NICE,
			NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);

	// Add a check fro n_rep to be >0 
	if (n_rep < 1 ) {gmx_fatal(FARGS,"The number of replicas should be larger than zero");}

	// As we need to generate and ensemble 
	// we need to output multiple xtcs, right now not implemented
	//out_file=opt2fn("-o",NFILE,fnm);
	//fprintf(stderr,"Will write %s \n",out_file);
	//trxout = open_xtc(out_file,"w");

	// Read RDC file
	rdcfile = opt2fn("-r",NFILE,fnm);
	if(rdcfile) {f_rdc = ffopen(rdcfile,"r");}
	while(fgets(line, 255, f_rdc) != NULL)
	{
		n_rdc++;
	}

	if (bVerbose) {
		printf("There are %d RDCs from %s\n",n_rdc,rdcfile);
	}
	rewind(f_rdc);

	// Allocate data for rdcs
	if((rdc = (rdc_str_e *) calloc(n_rdc, sizeof(rdc_str_e))) == NULL) {
		printf("\nERROR, cannot allocate rdc_str\n"); exit(1);
	}
	for (i=0; i<n_rdc; i++) {
		snew(rdc[i].comp,n_rep);
	}
	// Read in the RDCs 
	// Crude parsing, no comments or extra lines allowed
	j=0;
	while(fgets(line, 255, f_rdc) != NULL)
	{
		sscanf(line, "%d %d %f %f", &rdc[j].ind[0], &rdc[j].ind[1], &rdc[j].exp, &rdc[j].gyr);
		j++;
	}

	read_tps_conf(ftp2fn(efTPS, NFILE, fnm), buf, &top, &ePBC, &xp, NULL, box, TRUE);
	natoms = read_first_x(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &t, &x, box);

	do {

		// Find best orientation for the frame, compute RDCs and store them for averaging
		find_best_orient(rdc, n_rep, n_rdc, natoms, x, max_iter, mc_half_step, inv_T, bVerbose ) ; 
		//write_xtc(trxout, natoms, step, t, box, x, 1000);
		nframes++;

	} while (read_next_x(oenv, status, &t, natoms, x, box));

	//close_xtc(trxout);

	// Perform the ensemble average of RDCs and the Q-factor of them
	av_qfact=compute_av_qfact( rdc, n_rdc, nframes ); 
	printf("Trajectory averaged Q-factor %f \n", av_qfact);

	// Output computed vs experimental RDCs
	backfile = opt2fn("-q",NFILE,fnm);
	if(backfile) {
		f_back = ffopen(backfile,"w");
		write_backrdc(f_back, rdc, n_rdc, av_qfact);
	}

	return 0;

}

