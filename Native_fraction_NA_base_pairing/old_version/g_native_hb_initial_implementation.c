//
// Compute HB between atoms given in a file
// that needs to be supplied to the program
//
/* This line is only for CVS version info */
static char *SRCID_template_c = "$Id: g_hb_native, 0.1  2009/01/21 12:15 Portella $";

#include "statutil.h"
#include <string.h>
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "copyrite.h"
#include "statutil.h"
#include "tpxio.h"
#include "index.h"
#include "gmxfio.h"
#include <stdio.h>
#include <futil.h>
#include "rmpbc.h"
#include "pbc.h"
#include "mshift.h"
#include "do_fit.h"
#include "pdbio.h"
#include "mdrun.h"
#include "xvgr.h"


#define PI      3.14159265358979323846   /* PI = 4*atan(1)       */
#define BOLTZ   0.00831451               /* kJ / (mol*K)         */


t_atoms  *top_atoms=NULL;


FILE *in_p;
FILE *out;
FILE *fn;

real compute_native_hb(int *data[], int counter, rvec *x)
{

	int i,j,k;
	real native_hb=0;
	real *fractions;
	snew(fractions, counter);
    double nn=8, mm=14;
    rvec rij_d,rij, sij;
    real mod_rij_d;
    real mod_rij, mod_sij, rdist, r_alpha;
    double nhbond, r6dist, r12dist, num, iden;
    double r6ang, r12ang;
    double rijfmod;
    double cutang=120. * (2.0*3.1415927/360.);
    double d_contrib=0;
    double a_contrib=0;
    double alpha=0;
    double t1, t2, t3, t4, t5, t6;
    double caa, cbb, cab, ccc;
    double ac;


	for(i=0; i< counter; i++){
		fractions[i]=0;
	}

	for(i=0; i< counter; i++){

		rij[0] = x[data[1][i]-1][0]-x[data[2][i]-1][0];
		rij[1] = x[data[1][i]-1][1]-x[data[2][i]-1][1];
		rij[2] = x[data[1][i]-1][2]-x[data[2][i]-1][2];

		mod_rij_d  = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
		rdist= mod_rij_d/0.30;
		if(rdist>0.99998 && rdist<1.00002) rdist = 0.99999;       // to keep the function continuos
		r6dist = pow(rdist, nn);
		num = 1.-r6dist;                                          // numerator
		r12dist = pow(rdist, mm);
		iden = 1./(1.-r12dist);
		// num/idem is the distance contribution to hbond
		d_contrib =  num*iden;

		rij[0] = x[data[0][i]-1][0]-x[data[1][i]-1][0];
		rij[1] = x[data[0][i]-1][1]-x[data[1][i]-1][1];
		rij[2] = x[data[0][i]-1][2]-x[data[1][i]-1][2];

		mod_rij  = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]) ;

		sij[0] = x[data[1][i]-1][0]-x[data[2][i]-1][0];
		sij[1] = x[data[1][i]-1][1]-x[data[2][i]-1][1];
		sij[2] = x[data[1][i]-1][2]-x[data[2][i]-1][2];

		mod_sij = sqrt(sij[0]*sij[0]+sij[1]*sij[1]+sij[2]*sij[2]);
		t1 = rij[0];
		t2 = rij[1];
		t3 = rij[2];
		t4 = sij[0];
		t5 = sij[1];
		t6 = sij[2];
		caa = t1*t1+t2*t2+t3*t3;
		cab = t1*t4+t2*t5+t3*t6;
		cbb = t4*t4+t5*t5+t6*t6;
		ccc = gmx_invsqrt(caa*cbb);

		ac = cab*ccc; // cosinus de teta
		alpha = acos(ac);
		r_alpha= alpha/0.52;
		if(r_alpha>0.99998 && r_alpha<1.00002) r_alpha = 0.99999;       // to keep the function continuos
		r6ang = pow(r_alpha, nn);
		num = 1.-r6ang;                                          // numerator
		r12ang = pow(r_alpha, mm);
		iden = 1./(1.-r12ang);
		a_contrib= num*iden;

		fractions[i] = d_contrib*a_contrib;


	}
	
	int ib4=0;
	float vb4=0;
	real av=0;
	k=0;
	for(i=0; i<counter; i++){

		if(i==0){

			ib4=data[3][i];
			t1=fractions[0];

		}else if (i==(counter-1)){
			t1+=fractions[i];
			av = t1/(int)data[4][i];
			native_hb += av;
		}else{

			if(data[3][i]==ib4){

				t1+=fractions[i];
				ib4=data[3][i];

			}else{

				av = t1/(int)data[4][i-1];
				native_hb += av;
				t1=0;
				ib4=data[3][i];
				t1=fractions[i];

			}	
		}

	}

	sfree(fractions);
	return native_hb;

}
 


real dist3dpbc(real *point, real *cent, real x_box, real y_box, real z_box)
{

        real dist;
        real A=0;
        real B=0;
        real C=0;
        real dis_x=0;
        real dis_y=0;
        real dis_z=0;

        dis_x = cent[0] - point[0];
        dis_y = cent[1] - point[1];
        dis_z = cent[2] - point[2];

	// Aply PBC
        while  (dis_x > x_box*0.5){ dis_x= dis_x-x_box;}
        while  (dis_x < x_box*(-0.5)){ dis_x= dis_x+x_box;}
        while  (dis_y > y_box*0.5){ dis_y= dis_y-y_box;}
        while  (dis_y < y_box*(-0.5)){ dis_y= dis_y+y_box;}
        while  (dis_z > z_box*0.5){ dis_z= dis_z-z_box;}
        while  (dis_z < z_box*(-0.5)){ dis_z= dis_z+z_box;}

        dist = (sqrt( dis_x*dis_x + dis_y*dis_y + dis_z*dis_z )) ;

        return(dist);
}

real sqrt_dist3dpbc(real *point, real *cent, real x_box, real y_box, real z_box)
{
        real sqrt_dist;
        real A=0;
        real B=0;
        real C=0;
        real dis_x=0;
        real dis_y=0;
        real dis_z=0;

        dis_x = cent[0] - point[0];
        dis_y = cent[1] - point[1];
        dis_z = cent[2] - point[2];

        while  (dis_x > x_box*0.5){ dis_x= dis_x-x_box;}
        while  (dis_x < x_box*(-0.5)){ dis_x= dis_x+x_box;}
        while  (dis_y > y_box*0.5){ dis_y= dis_y-y_box;}
        while  (dis_y < y_box*(-0.5)){ dis_y= dis_y+y_box;}
        while  (dis_z > z_box*0.5){ dis_z= dis_z-z_box;}
        while  (dis_z < z_box*(-0.5)){ dis_z= dis_z+z_box;}

        sqrt_dist = ( dis_x*dis_x + dis_y*dis_y + dis_z*dis_z ) ;

        return(sqrt_dist);
}



int main(int argc,char *argv[])
{

	matrix rotation_mat;
	matrix inv_rotation_mat;
	t_topology top;
	char line[255];
	int        ePBC;
	char       title[STRLEN];
	t_trxframe fr;
	rvec       *xtop;
	matrix     box;
	int        status;
	int        flags = TRX_READ_X;
	static int n=1;
	bool    bIndex;
	bool	bTPS;
	int 	i,j,k,l,m,r,d, a;
	int     iselect, ifit;
	atom_id *ind_select,*ind_fit;
	char    *gn_select,*gn_fit;
	bool    bFit=FALSE;
    real    *w_rls=NULL;
	t_pbc   *pbc;
	int 	counter=0;
	bool bDebug=FALSE;
	int  *data[5];
	
	bIndex=TRUE;


	//We should write out the density such that it can be re-analyzed afterwards...

	static char *desc[] = {
		"This g_program computes the number of native pairs hb. [BR]",
		"[BR]",
		"We use a smooth cut-off criteria based on the fulfillment of three criteria:[BR]",
		"[BR]",
		"------------------------------------------------------------------------------[BR]",
		" IMPORTANT!: So far PBC are not taken into account, so feed me well![BR]",
		" IMPORTANT!: So far RNA is not taken into account, add them to hardcoded lib ![BR]",
		"------------------------------------------------------------------------------[BR]",
		"[BR]",
		"	",
		""
	};

    t_pargs pa[] = {
        {"-fit", TRUE, etBOOL, {&bFit}, "Fit the selection."},
    };


	t_filenm fnm[] = {
		{ efTPS,  NULL,  NULL, ffREAD },   /* this is for the topology */
		{ efTRX, "-f", NULL, ffREAD },      /* and this for the trajectory */
		{ efXVG, "-o",  "native_hb",   ffWRITE },  /* output the info, time vs stacking */
		{ efDAT, "-p" , NULL, ffREAD },      /* and this for the index file*/
	};

#define NFILE asize(fnm)

	CopyRight(stderr,argv[0]);

	parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,
			NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);


	//bIndex=ftp2bSet(efNDX,NFILE,fnm);
	bTPS = ftp2bSet(efTPS,NFILE,fnm) ;

	if ( bTPS ) {
		read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&xtop,NULL,box,TRUE);
		//sfree(xtop);
		top_atoms=&top.atoms;

	} else {
		printf("\nSTOP: I need some sort of topology or pdb structure file\n\n");
		exit(1);
	}

	if (bIndex && bFit) {
		printf("Select group for the analysis:\n");
		get_index(top_atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&iselect,&ind_select,&gn_select);
		printf("Now select group for the fitting\n");
		get_index(top_atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&ifit,&ind_fit,&gn_fit);
		snew(w_rls,top.atoms.nr);
		for(i=0; (i<ifit); i++) { w_rls[ind_fit[i]]=top_atoms->atom[ind_fit[i]].m; } // this are the atoms to fit 
	}

	in_p=ffopen(opt2fn("-p",NFILE,fnm),"r");

    while(fgets(line, 255, in_p) != NULL)
	{
		counter++;
	}

    if((data[0] = ( int *) calloc(counter, sizeof( int))) == NULL) {printf("\nERROR, cannot allocate data\n"); exit(1);}
    if((data[1] = ( int *) calloc(counter, sizeof( int))) == NULL) {printf("\nERROR, cannot allocate data\n"); exit(1);}
    if((data[2] = ( int *) calloc(counter, sizeof( int))) == NULL) {printf("\nERROR, cannot allocate data\n"); exit(1);}
    if((data[3] = ( int *) calloc(counter, sizeof( int))) == NULL) {printf("\nERROR, cannot allocate data\n"); exit(1);}
    if((data[4] = ( int *) calloc(counter, sizeof( int))) == NULL) {printf("\nERROR, cannot allocate data\n"); exit(1);}
	// Seems 0->2 are pointers to donor-h-acceptor
	// the last two... one might be residue number and the other a weight factor??

	rewind(in_p);
	j=0;
    while(fgets(line, 255, in_p) != NULL)
	{
		 sscanf(line, "%d %d %d %d %d", &data[0][j],&data[1][j], &data[2][j],&data[3][j],&data[4][j]);
		 j++;
	}

	for(i=0; i<counter;i++){
		 fprintf(stdout, "%d %d %d %d %d\n", data[0][i],data[1][i], data[2][i],data[3][i], data[4][i]);
	
	}
	

	printf("\nFound %d values\n", counter);

	 if (ePBC != epbcNONE) {
		snew(pbc,1);
	 }else {
		 pbc = NULL;	
	 }

	read_first_frame(&status,ftp2fn(efTRX,NFILE,fnm),&fr,flags);
	fflush(stdout);
	printf("\nFound this box size: %f %f %f\n",fr.box[0][0],fr.box[1][1], fr.box[2][2]);


	if (pbc) {
		set_pbc(pbc,ePBC,box);
		rm_pbc(&top.idef,ePBC,top.atoms.nr,box,fr.x,fr.x);
	}
    rvec x_shift;
	rvec xcm;
	
	out=ffopen(opt2fn("-o",NFILE,fnm),"w");

	real native_hb=0;
	do {

		native_hb=compute_native_hb(data, counter, fr.x);
		fprintf(out,"%f %f\n", fr.time, native_hb);

	} while(read_next_frame(status,&fr));

	thanx(stderr);
	ffclose(out);
	return 0;

}

