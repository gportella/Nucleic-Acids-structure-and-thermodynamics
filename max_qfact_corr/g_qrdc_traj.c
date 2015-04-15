// GPC 18/03/2015 g_qrdc_traj v 0.999
// Code to back-compute RDCs the theta-method way. Finds the optimal correlation
// for computed vs exp RDCs by x/y rotations, then back-calculates the RDCs in
// the best orientation. Based on RDC colvar from plumed 
//
// Easily fixable caveats : No pbc taken into account, whole molecules only
// "Feature" : Only one scaling for all type of bonds
//
// TODO: 
// -- Better parsing and checking of the input RDC file
// -- Add more error checks for routines that could fail
//
// RDC input file free format, but no comments or extra lines allowed. The
// format is 
//
// atom_index atom_index RDC gyromagnetic_ratio
// 
// Weights format is simply the column with the weight of the frame
// its row number equals the number of frame (no matter time resolution)
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
#include <omp.h>

#define Const (0.3356806)
#define PI      3.14159265358979323846

typedef struct{
	int ind[2];		// index for each atom in the bond / atom pair
	real exp;		// experimental RDC
	real *comp;		// computed RDC
	rvec *pos[2];	// xyz coords of atoms in atom pair
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

void read_traj(rdc_str_e *rdc, int n_rdc,rvec *x, int nfr ) 
{

	int i,j;
	int a1, a2;
	for (i=0; i<n_rdc; i++){
		a1 = rdc[i].ind[0]-1 ;
		a2 = rdc[i].ind[1]-1 ;
		copy_rvec(x[a1],rdc[i].pos[0][nfr]);
		copy_rvec(x[a2],rdc[i].pos[1][nfr]);
	}

}

void shuffle_rdc( rdc_str_e *rdc, int n_rdc)
{
	//Fisher-Yates 
	int n = n_rdc;  
	int k=0; 
	int i=0;
	real x_sup;

	while (n > 1) {  
		n--;  
		k = randint(n) ;
		x_sup = rdc[k].exp ;
		rdc[k].exp = rdc[n].exp ;
		rdc[n].exp = x_sup ;
	}

}

void shuffle_traj(rdc_str_e *rdc, real *we, int  n_rdc, int nfr) 
{
	//Fisher-Yates 
	int n = nfr;  
	int k=0; 
	int i=0;
	real sup;
	rvec x_sup; 

	while (n > 1) {  
		n--;  
		k = randint(n) ;
		for(i=0; i< n_rdc; i++){
			//exchange w random pair
			copy_rvec(rdc[i].pos[0][k], x_sup);
			copy_rvec(rdc[i].pos[0][n], rdc[i].pos[0][k]);
			copy_rvec(x_sup, rdc[i].pos[0][n]);

			copy_rvec(rdc[i].pos[1][k], x_sup);
			copy_rvec(rdc[i].pos[1][n], rdc[i].pos[1][k]);
			copy_rvec(x_sup, rdc[i].pos[1][n]);

			sup = we[k];
			we[k] = we[n];
			we[n] = sup;
		}
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



real compute_av_qfact( rdc_str_e *rdc, int n_rdc, int wp )
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
		rdc[i].av_comp /= wp;
		sum+=(rdc[i].exp-rdc[i].av_comp)*(rdc[i].exp-rdc[i].av_comp);
		Q+=(rdc[i].exp*rdc[i].exp);
	}
	q_fact = sqrt(sum/Q); 
	return q_fact;
}

real compute_rdc_qfact(rdc_str_e *rdc, int n_rep, int n_rdc, int st, int n_len, 
		rvec **x, real reg_mes, real *lw,  bool bCorrel)
{
	//
	// Returns either the correlation between computed and experimental RDC or the
	// Q-factor between them. Recomended use is with scaling factors of 1 for
	// bCorrel == True, and with the actual scaling factor otherwise.  
	// Addapted from plumed's colvar RDC.
	//
	int a1=0;
	int a2=0;
	int i,j,r;
	int ct=0;
	int k=0;
	real scx=0;
	real scx2=0;
	real scy=0;
	real scy2=0;
	real scxy=0;
	real sum=0;
	real Q=0;
	real S=0;
	real corr_qfact=0;
	real r_norm=0;
	real xx,yy,zz,d ;

	ct = 2 * n_rdc;
	r_norm = 1. / (n_len);
	real norm_weight = 1. / (real) n_rep;

	real **weight;
	// fixed grid for the entropy penalty based on 
	// orientation of an arbitray RDC, which we select
	// to be the last one in the list. It is not a measure
	// of the number conformation, but rotations. We would
	// like to minimize having too many different rotations
	// to explain the same set of RDCs. It should lower the 
	// score for large number of replicas.
	int grid_theta=18;
	int grid_phi=36;
	real dtheta=180./(real)grid_theta;
	real dphi=360./(real)grid_phi;
	snew(weight,grid_theta);
	for (i=0; i< grid_theta; i++) { snew(weight[i],grid_phi); }

	//we set the ensemble averaged RDC to zero 
	for (i=0; i< n_rdc; i++) { rdc[i].r_comp = 0; }

	k=0;
	int chunk=10;
	for (r=0; r< n_len; r++){
		int k=0;
		// rdc vectors are stored sequentially 
		// ever second is another vector
		for (i=0; i < ct ; i = i+2){
			// this could be changed to account for PBCs 
			xx = x[r][i][XX] - x[r][i+1][XX] ;
			yy = x[r][i][YY] - x[r][i+1][YY] ;
			zz = x[r][i][ZZ] - x[r][i+1][ZZ] ;
			d = sqrt( xx*xx + yy*yy + zz*zz );
			real d2   = d*d;
			real d3   = d2*d;
			real id3  = 1./d3;
			real max  = -Const*rdc[k].scale*rdc[k].gyr;
			real dmax = id3*max;
			real cos_theta = zz / d ;
		//	rdc[k].comp[r] = 0.5*dmax*(3.*cos_theta*cos_theta-1.);
			// replica averaged RDC, here is where we include the weight
			rdc[k].r_comp += r_norm*lw[r]*0.5*dmax*(3.*cos_theta*cos_theta-1.) ;
			k++;
		}
		if(reg_mes != 0 ){
			// we compute the hist of orientations of the last vector
			// domain [0,pi], both included
			real theta = acos(zz/d);
			// we shift the domain of atan2 to [0,2pi]
			// the angle of course will be wrong, but 
			// we don't care, we don't want to use ifs 
			real phi = atan2(yy,xx)+PI;
			// the last term should avoid problems 
			// if theta == pi or phi == 2pi 
			a1 = (int) floor(theta/(dtheta*DEG2RAD)) - (int) floor(theta/PI); 
			a2 = (int) floor((phi)/(dphi*DEG2RAD)) - (int) floor(phi/(2*PI)); 
			//populate the grid
			weight[a1][a2] += norm_weight; 
		}
	}
	if(reg_mes != 0 ){
		// we collect the entropy penalty
		for(i=0; i<grid_theta-2; i++){
			for(j=0; j<grid_phi-2; j++){
				S += weight[i][j]*log(n_rep*weight[i][j]+1);
			}
		}
		// scale by given factor reg_mes
		S *= reg_mes;
	}
	// statistics for correlation coefficient
	for (i=0; i< n_rdc ; i++){
		scx  += rdc[i].r_comp;
		scx2  += rdc[i].r_comp*rdc[i].r_comp;
		scy  += rdc[i].exp;
		scy2  += rdc[i].exp*rdc[i].exp;
		scxy += rdc[i].exp*rdc[i].r_comp;
	}
	real num = (n_rdc*n_len)*scxy - scx*scy;
	real idevx = 1./sqrt((n_rdc*n_len)*scx2-scx*scx);
	real idevy = 1./sqrt((n_rdc*n_len)*scy2-scy*scy);
	// the correlation
	corr_qfact = num * idevx * idevy;

	if (!bCorrel)   {
		// Compute q-factor
		// We don't want to mix up the penalty for the MES
		S = 0;
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

	for (i=0; i< grid_theta; i++) { sfree(weight[i]); }
	sfree(weight);

	// slow...  How can i get rid of the if?
	// corr_qfact is not always positive
	if( corr_qfact > 0 ) { 
		return (corr_qfact-S);
	}else{
		return (corr_qfact);
	}
}

void find_best_orient(rdc_str_e *rdc, int n_rep, int n_rdc, int st, int end, int max_iter, 
		real mc_half_step, real inv_T, real anneal_w, real reg_mes, real *we, bool bVerbose)
{
	//
	// Find the alignment that best correlates the computed RDCs and the
	// experimental ones. It uses a Metropolis Monte Carlo search. Parameters are
	// tunable, but defaults should be okeish. After the alignment is found,
	// computes the slope of the correlation (aka scaling) and computes the
	// Q-factor. The rdc_str_e structure containts the best rdcs per frame, and
	// a cumulative sum, such that we can perform the average later on. 
	//
	int i, j ,k, r, ct, acc, count, xrot;
	int rest=1;
	int annealing = 0;
	int r_rot=0;
	int n_len=0;
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
	real *lw;
	rvec **xens;
	rvec *xnew;
	rvec **xnew_e;
	rvec **xopt_e;
	bool bCorrel=TRUE;

	// this should be the same as n_rep in this implemention
	// where we do not have chuncks of diff size (ends dropped)
	n_len = end - st;
	//each rdc has 2 sets of coordinates
	ct = n_rdc*2;
	// proposed rotation of one frame
	snew(xnew,ct);
	// optimized coords all frames
	snew(xopt_e,n_rep);
	// proposed MC step all frames
	snew(xnew_e,n_rep);
	// accepted MC step all frames
	snew(xens,n_rep);
	// local weights
	snew(lw,n_rep);

	// create initial replicas of the frame
	for(r=0; r<n_rep; r++){
		snew(xens[r],ct);	
		snew(xnew_e[r],ct);	
		snew(xopt_e[r],ct);	
		for(j=0; j < n_rdc; j++){
			copy_rvec(rdc[j].pos[0][st+r],xens[r][k]);
			copy_rvec(rdc[j].pos[1][st+r],xens[r][k+1]);
			copy_rvec(rdc[j].pos[0][st+r],xnew_e[r][k]);
			copy_rvec(rdc[j].pos[1][st+r],xnew_e[r][k+1]);
			copy_rvec(rdc[j].pos[0][st+r],xopt_e[r][k]);
			copy_rvec(rdc[j].pos[1][st+r],xopt_e[r][k+1]);
			// weight
			lw[r]=we[st+r];

			k=k+2;
		}
		k=0;
	}

	// set the scaling to 1 for correlation
	for (i=0; i< n_rdc; i++) { rdc[i].scale = 1  ; }
	// initiate qfact
	corr_b4 = compute_rdc_qfact(rdc, n_rep, n_rdc, st, n_len, xens, reg_mes, lw, bCorrel);
	max_corr = corr_b4;
	real orig_inv_T = inv_T;
	count=0;
	acc=0;
	do {
		count++;
		// do annealing to optimize the acceptance, such that large rotations
		// are favoured at the beginning but not at the end, presumably near the max
		annealing++;
		inv_T = orig_inv_T + log(anneal_w*annealing*orig_inv_T+1);
		//select which one you want to rotate 
		r_rot = randint(n_rep) ;
		// generate random rotation [M,N] range
		// M + rand() / (RAND_MAX / (N - M + 1) + 1)
		xr = -mc_half_step + ((real) rand() / ( (real) RAND_MAX / (mc_half_step*2)+1 )) ;
		yr = -mc_half_step + ((real) rand() / ( (real) RAND_MAX / (mc_half_step*2)+1 )) ;
		rotx = DEG2RAD * xr;
		roty = DEG2RAD * yr;
		// propose new rotation
		rotate(ct, xens[r_rot], xnew, rotx, roty, 0);
		//copy xnew to its place in xnew_e
		for(i=0; i<ct; i++){ copy_rvec(xnew[i],xnew_e[r_rot][i]); }
		// Try qfact rotated ala Metropolis
		corr = compute_rdc_qfact(rdc, n_rep, n_rdc, st, n_len, xnew_e, reg_mes, lw, bCorrel);
		ran_met = ((real) rand() / (real) RAND_MAX )*(1);
		if(bVerbose){
			printf("Ef. Corr: %5.3f, Max Ef. corr: %5.3f, acc. rat. %5.3f, c %8d inv_T %8.6f\r", corr, 
					max_corr, (real) acc / (real) count, count, inv_T);
			fflush(stdout);
		}

		if (exp(-(corr_b4-corr)/inv_T) < ran_met  || corr > corr_b4){
			if (corr > max_corr){
				max_corr = corr ; 
				// We keep it if it improves
				for(r=0; r< n_rep ; r++){
					for(i=0;i<ct; i++){
						copy_rvec(xnew_e[r][i],xopt_e[r][i]);
					}
				}
			}
			// Copy it for next iteration
			for(i=0;i<ct; i++){copy_rvec(xnew_e[r_rot][i],xens[r_rot][i]);}
			corr_b4 = corr ; 
			acc++;
		}
	} while ( max_corr < 0.999  && count < max_iter );

	// rdc structure has been updated with newest computed values
	corr=compute_rdc_qfact(rdc, n_rep, n_rdc, st, n_len, xopt_e, reg_mes, lw, bCorrel);
	if (bVerbose) {
		printf ("\nMax. correlation after MMC is %f (%d iterations)\n", corr, count);
	} 
	k=0;
	//copy opt to rdc structure
	for(r=0; r<n_rep; r++){
		for(j=0; j < n_rdc; j=j+2){
			copy_rvec(xopt_e[r][k],rdc[j].pos[0][st+r]);
			copy_rvec(xopt_e[r][k+1],rdc[j].pos[1][st+r]);
			k=k+2;
		}
		k=0;
	}
	// fill up arrays for fitting
	snew(fx,n_rdc);
	snew(fy,n_rdc);
	for(i=0; i<n_rdc;i++){
		fx[i]=rdc[i].r_comp;
		fy[i]=rdc[i].exp;
	}
	lsq_y_ax( n_rdc, fx, fy, &slope); // fitting  y=ax
	//printf("Slope if %f\n", slope);
	sfree(fx);
	sfree(fy);
	// the slope is the scaling factor
	for (i=0; i< n_rdc; i++) { rdc[i].scale = slope ; }
	// Now we do want to compute a Q-factor
	// The rdc addition for average calculation takes place
	bCorrel=FALSE;
	qfact=compute_rdc_qfact(rdc, n_rep, n_rdc, st, n_len, xopt_e, reg_mes, lw, bCorrel);
	if (bVerbose){
		printf("Per frame Q-fact is %f \n",qfact);
	}

	// free everything local
	for(r=0; r< n_rep; r++){
		sfree(xnew_e[r]);
		sfree(xopt_e[r]);
		sfree(xens[r]);	
	}
	sfree(xnew);
	sfree(xopt_e);
	sfree(xnew_e);
	sfree(xens);
	sfree(lw);

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
		"between the ensemble of structures and the experimental values.[BR]"
		"The user should select an appropriate number of frames to be.[BR]"
		"time averaged. If this number is too large, it leads to overfitting.[BR]"
		"[BR]"
		"RDC input file free format, but no comments or extra lines allowed.[BR]"
		"the format is [BR] " 
		"[BR]"
		"atom_index atom_index RDC gyromagnetic_ratio[BR]"
		"[BR]"
		"Weights format is simply the column with the weight of the frame[BR]"
		"its row number equals the number of frame (no matter time resolution)[BR]"
		"We only use weights if given by the user."
	};

	FILE				*f_rdc;
	FILE				*f_back;
	FILE				*f_w;
	static int	max_iter=50000;
	static int	n_rep=28;
	static real	mc_half_step=5;
	static real	inv_T=0.0001;
	static real	reg_mes=0.1;
	static real	anneal_w=0.001;
	static bool	bVerbose=FALSE;
	static bool	bRandom=FALSE;

	/* Command-line arguments */
	t_pargs          pa[] = {
		{"-mi", FALSE, etINT, {&max_iter},
			"Maximum number of iterations MC search"},
		{"-ms", FALSE, etREAL, {&mc_half_step},
			"MC half step size, in degrees"},
		{"-mt", FALSE, etREAL, {&inv_T},
			"MC weight factor / inverse temperature"},
		{"-nr", FALSE, etINT, {&n_rep},
			"Number of replicas"},
		{"-ma", FALSE, etREAL, {&anneal_w},
			"Scaling factor s for annealing, inv_T(t) = log(s*time*inv_T(0)+1)"},
		{"-mes", FALSE, etREAL, {&reg_mes},
			"Scaling factor s for Max. Entropy regularisation s*Sum(w*log(wN+1))"},
		{"-ran", TRUE, etBOOL, {&bRandom},
			"Randomize the RDCs, to do validation."},
		{"-v", TRUE, etBOOL, {&bVerbose},
			"Be verbose, output correlations for each replica to std output."},
	};

	/* Output files */
	t_filenm            fnm[] = {
		{ efTRX, "-f", NULL, ffREAD },
		{ efTPS, NULL, NULL, ffREAD },
		{ efDAT, "-r", "rdc.dat", ffREAD },
		{ efDAT, "-w", "weights.dat", ffREAD },
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
	real			*we; //weights
	int				natoms,i,j,kkk;
	int				ePBC;
	int				nframes=0;
	int				n_rdc=0;
	int				step=1;
	int				st=0;
	int				wp=0;
	int				endf=0;
	int				n_w=0;
	char			buf[256];
	const char		*rdcfile=NULL;
	const char		*backfile=NULL;
	const char		*out_file=NULL; 
	const char		*wfile=NULL; 
	char			line[255];
	bool			bWeights=FALSE;


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

	if (n_rep > n_rdc  ) {gmx_warning("The number of replicas should be smaller than"
			" the number of RDCs");}

	// Allocate data for rdcs
	int rdc_alloc=10000;
	int rdc_alloc_orig= rdc_alloc;
	snew(rdc,n_rdc);
	for (i=0; i<n_rdc; i++) { 
		snew(rdc[i].comp,rdc_alloc);
		snew(rdc[i].pos[0],rdc_alloc);
		snew(rdc[i].pos[1],rdc_alloc);
	}

	// Read in the RDCs 
	// Crude parsing, no comments or extra lines allowed
	j=0;
	while(fgets(line, 255, f_rdc) != NULL)
	{
		sscanf(line, "%d %d %f %f", &rdc[j].ind[0], &rdc[j].ind[1], &rdc[j].exp, &rdc[j].gyr);
		j++;
	}

	// Read atoms from pdb file
	read_tps_conf(ftp2fn(efTPS, NFILE, fnm), buf, &top, &ePBC, &xp, NULL, box, TRUE);
	natoms = read_first_x(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &t, &x, box);


	do {
		// Read the trajectory
		// We don't want to read the traj more than once, therefore we 
		// allocate memory in advanced, and we change its size if needed
		if (nframes >= rdc_alloc) { 
			rdc_alloc +=  rdc_alloc_orig ; 
			for (i=0; i<n_rdc; i++) { 
				srenew(rdc[i].comp,rdc_alloc); 
				srenew(rdc[i].pos[0],rdc_alloc); 
				srenew(rdc[i].pos[1],rdc_alloc); 
			}
		}
		read_traj(rdc, n_rdc, x, nframes) ; 
		nframes++;

	} while (read_next_x(oenv, status, &t, natoms, x, box));
	printf("Read %d frames\n", nframes);

	// Read weight if needed
	bWeights=opt2bSet("-w",NFILE,fnm);
	if (bWeights){
		wfile = opt2fn("-w",NFILE,fnm);
		if(wfile) {f_w = ffopen(wfile,"r");}
		while(fgets(line, 255, f_w) != NULL) { n_w++; } 
		rewind(f_w);
		if( n_w != nframes ) {
			printf("Number of frames %d , number of weights %d\n", nframes, n_w);
			gmx_fatal(FARGS,"Sorry, number of frames and number of weights does not match"); 
		}
	}	
	j=0;
	// allocate for weights, read them in
	// or else make them equal 
	snew(we,nframes);
	if (bWeights){
		while(fgets(line, 255, f_w) != NULL)
		{
			sscanf(line, "%f", &we[j]);
			j++;
		}
	}else{
		for(i=0; i<nframes; i++){ we[i] = 1.0 / nframes ; }
	}

	//initiate PRNG
	srand(time(NULL));
	// now shuffle frames in trajectory, 
	// the weights reorder accordingly
	shuffle_traj(rdc, we, n_rdc, nframes) ; 
	if (bRandom){
		shuffle_rdc(rdc, n_rdc) ; 
	}

	// Divide the trajectory in chunks
	endf=n_rep;
	wp = floor((real)nframes / (real) n_rep);

	// we throw away what's left after ceil, like this
	for(i=0; i<wp; i++){
		st = i*n_rep;
		endf = st+n_rep;
		find_best_orient(rdc, n_rep, n_rdc, st, endf, max_iter, 
				mc_half_step, inv_T, anneal_w, reg_mes, we, bVerbose ) ; 
		if(!bVerbose){
			printf("\rDone %5.2f %%", 100* ( (real) i / (real) wp));
			fflush(stdout);
		}
	}

	// Perform the ensemble average of RDCs and the Q-factor of them
	av_qfact=compute_av_qfact( rdc, n_rdc, wp ); 
	printf("\nTrajectory averaged Q-factor %f \n", av_qfact);

	// Output computed vs experimental RDCs
	backfile = opt2fn("-q",NFILE,fnm);
	if(backfile) {
		f_back = ffopen(backfile,"w");
		write_backrdc(f_back, rdc, n_rdc, av_qfact);
	}

	return 0;

}

