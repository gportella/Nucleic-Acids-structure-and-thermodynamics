//
// Compute HB between atoms given in a file
// that needs to be supplied to the program
// 


/* This line is only for CVS version info */
static char *SRCID_template_c = "$Id: g_stacking, 0.1  2009/01/21 12:15 Portella $";

#include <gromacs/copyrite.h>
#include <gromacs/filenm.h>
#include <gromacs/macros.h>
#include <gromacs/pbc.h>
#include <gromacs/smalloc.h>
#include <gromacs/statutil.h>
#include <gromacs/rmpbc.h>
#include <gromacs/xvgr.h>
#include <gromacs/gmx_fatal.h>

#include <gromacs/mdrun.h>
#include <gromacs/force.h>
#include <gromacs/nonbonded.h>
#include <gromacs/tpxio.h>
#include <gromacs/mshift.h>
#include <gromacs/vec.h>
#include <gromacs/nrnb.h>
#include <gromacs/typedefs.h>
#include <string.h>
#include <stdbool.h>

#define PI      3.14159265358979323846   /* PI = 4*atan(1)       */

typedef struct{
	int resnum;
	int *don;
	int ndon;
	int *h;
	int nh;
	int *acc;
	int nacc;
} base;

t_atoms  *top_atoms=NULL;
char *pur[]={ "N1","C2","N3","C4","C5","C6","N7","C8","N9"};
char *pyr[]={ "N1","C2","N3","C4","C5","C6"};

char *bpur[]={"DG5", "DG", "DG3","DA", "DA3", "DA5", "RG5", "RG", "RG3","RA", "RA3", "RA5", "A","G", "AMO", "GMO", "A5E", "A3E"  };
char *bpyr[]={"DC", "DC3", "DC5", "DT", "DT3", "DT5","RT5", "RC", "RC3", "RC5", "RT", "RT3", "RT5", "C","T", "U", "TMO", "CMO", "PMO", "DP", "DP5", "DP3", "T5E" };
char *na[]={ "DG5", "DG", "DG3","DA", "DA3", "DA5", "RG5", "RG", "RG3","RA", "RA3", "RA5", "A","G", "AMO", "GMO  ", "A5E", "A3E", "DC", "DC3", "DC5", "DT", "DT3", "DT5","RT5", "RC", "RC3", "RC5", "RT", "RT3", "RT5", "C","T", "U"  , "TMO", "CMO", "PMO", "DP", "DP5", "DP3", "T5E"};

char *bA[]={ "N1","C2","N3","C4","C5","C6","N7","C8","N9"};
char *bG[]={ "N1","C2","N3","C4","C5","C6","N7","C8","N9"};
char *bT[]={ "N1","C2","N3","C4","C5","C6"};
char *bC[]={ "N1","C2","N3","C4","C5","C6"};

char *b_na[]= { "C2","C4","C5","C6","C7","C8","H1","H2","H21","H22","H3","H41","H42","H5","H6",
	        "H61","H62","H71","H72","H73","H8","N1","N2","N3","N4","N6","N7","N9","O2","O4","O6"};

char *cap_b[]={"DG5","DA5","DT5","DC5"};
char *cap_e[]={"DG3","DA3","DT3","DC3"};

char *hb_T_d[]={"N3"};
char *hb_T_h[]={"H3"};
char *hb_T_a[]={"O4"};

char *hb_A_d[]={"N6","N6"};
char *hb_A_h[]={"H61","H62"};
char *hb_A_a[]={"N1","N7"};

char *hb_C_d[]={"N4","N4"};
char *hb_C_h[]={"H41","H42"}; 
char *hb_C_a[]={"N3","O2"};

char *hb_G_d[]={"N2", "N2", "N1"}; 
char *hb_G_h[]={"H21", "H22","H1" }; 
char *hb_G_a[]={"N1", "N22"}; 

FILE *in_p;
FILE *out;
FILE *fn;

typedef enum { DA, DT, DC, DG } en_base;

int return_bp_type(char n1, char n2)
{
	// default is 0, which will mean no hb allowed, zero in the matrix
	// entry
	int ret=0;
	if ( ((n1=='a') && (n2=='t')) || 
		((n1=='t') && (n2=='a')) ){
		ret=1;
	}
	if ( ((n1=='c') && (n2=='g')) || 
		((n1=='g') && (n2=='c')) ){
		ret=2;
	}
	// if you want to take other configurations into account
	// add more conditions and enums in calling switch 

	return ret;

}

char return_case(char *name)
{
	char dd='f';
	if ( (strcmp(name,"DA")==0) || (strcmp(name,"DA3")==0) || 
			(strcmp(name,"DA5")==0) ||  (strcmp(name,"A5E")==0) || 
			(strcmp(name,"A3E")==0))  {
		dd='a';
	}
	if ( (strcmp(name,"DT")==0) || (strcmp(name,"DT3")==0) || 
		(strcmp(name,"DT5")==0) || (strcmp(name,"T5E")==0) )  {
		dd='t';
	}
	if ( (strcmp(name,"DC")==0) || (strcmp(name,"DC3")==0) || 
			(strcmp(name,"DC5")==0)  )  {
		dd='c';
	}
	if ( (strcmp(name,"DG")==0) || (strcmp(name,"DG3")==0) || 
			(strcmp(name,"DG5")==0) )  {
		dd='g';
	}

	return dd;
}

bool in_list(char *list[], int size, char *resnm )
{
	/* gives true if resnm occurs in the list of names provided */
	int i;

	for(i=0; (i<size); i++){
		if (strcmp(list[i],resnm) == 0){ return TRUE ;}
	}
	return FALSE;
}

void fill_base(t_topology top, base *b,int nres, char *d[],int nd, char *h[],int nh, char *a[], int na)
{
	int i,j;
	int dd=0;
	int hh=0;
	int aa=0;

	b->resnum = nres;
	for (j=0; j< nd; j++){
		for (i=0; i<top.atoms.nr; i++){
			if(	top.atoms.resinfo[top.atoms.atom[i].resind].nr == nres + 1 ){
				// comparing against donors
				if (strcmp(d[j], *(top.atoms.atomname[i])) == 0 )
				{
					b->don[dd] = i;
					dd++;
					if ( dd > b->ndon) {
						printf("Found more donors (%d) than expected (%d), %s %d\n",dd, b->ndon,  	
								*top.atoms.resinfo[nres].name, nres);
					}
				}

				// comparing against hydrogens
				if (strcmp(h[j], *(top.atoms.atomname[i])) == 0 )
				{
					b->h[hh] = i;
					hh++;
					if ( hh > b->nh) {
						printf("Found more donors (%d) than expected (%d), %s %d\n",hh, b->nh,  	
								*top.atoms.resinfo[nres].name, nres);
					}
				}
			}
		}
	}
	dd=0;
	hh=0;
	for (j=0; j< na; j++){
		for (i=0; i<top.atoms.nr; i++){
			if(	top.atoms.resinfo[top.atoms.atom[i].resind].nr == nres + 1 ){
				// comparing against acceptors
				if (strcmp(a[j], *(top.atoms.atomname[i])) == 0 )
				{
					b->acc[aa] = i;
					aa++;
					if ( aa > b->nacc) {
						printf("Found more donors (%d) than expected (%d), %s %d\n",aa, b->nacc,  	
								*top.atoms.resinfo[nres].name, nres);
					}
				}
			}
		}
	}
	aa=0;
}

void allocate_base(t_topology top, base *b,int nres, char *d[],int nd, char *h[],int nh, char *a[], int na)
{
	int i,j;
	b->ndon=0;
	b->nh=0;
	b->nacc=0;

	for (j=0; j< nd; j++){
		for (i=0; i<top.atoms.nr; i++){
			if(	top.atoms.resinfo[top.atoms.atom[i].resind].nr == nres + 1 ){
				if (strcmp(d[j], *(top.atoms.atomname[i])) == 0 )
				{
					b->ndon++;
				}
				// we should have as many d as h, so I use the same counter j
				if (strcmp(h[j], *(top.atoms.atomname[i])) == 0 )
				{
					b->nh++;
				}
			}
		}
	}
	for (j=0; j< na; j++){
		for (i=0; i<top.atoms.nr; i++){
			if(	top.atoms.resinfo[top.atoms.atom[i].resind].nr == nres + 1 ){
				if (strcmp(a[j], *(top.atoms.atomname[i])) == 0 )
				{
					b->nacc++;
				}
			}
		}
	}

	snew(b->don, b->ndon);
	snew(b->h, b->nh);
	snew(b->acc, b->nacc);
}

void populate_base(char *hb_T_d[],int s_T_d, char *hb_T_h[],int s_T_h, char *hb_T_a[], int s_T_a, char *hb_A_d[],int s_A_d, char *hb_A_h[],int s_A_h, char *hb_A_a[],int s_A_a, char *hb_C_d[],int s_C_d, char *hb_C_h[],int s_C_h, char *hb_C_a[], int s_C_a, char *hb_G_d[],int s_G_d,char *hb_G_h[],int s_G_h, char *hb_G_a[], int s_G_a, int natoms, t_topology top, base *base_bp, int **hb_max)
{

	int i, j, k, c_bases=0, sized_d=0;
	char *res1, *res2 ;
	char pp, p1, p2, bp_type;


	for (i=0; i<top.atoms.nres; i++)
	{
		res1=*top.atoms.resinfo[i].name;
		pp = return_case(res1);
		switch( pp )
		{
			case 'a':
				//printf ("Found DA, %s\n", resname);
				allocate_base(top, &base_bp[c_bases], i, hb_A_d, s_A_d, hb_A_h,s_A_h, hb_A_a, s_A_a);
				fill_base(top, &base_bp[c_bases], i, hb_A_d, s_A_d, hb_A_h, s_A_h, hb_A_a, s_A_a);
				c_bases++;
				break;
			case 't':
				//printf ("Found DT, %s\n", resname);
				allocate_base(top, &base_bp[c_bases], i, hb_T_d, s_T_d, hb_T_h,s_T_h, hb_T_a, s_T_a);
				fill_base(top, &base_bp[c_bases], i, hb_T_d, s_T_d, hb_T_h, s_T_h, hb_T_a, s_T_a);
				c_bases++;
				break;
			case 'c':
				//printf ("Found DC, %s\n", resname);
				allocate_base(top, &base_bp[c_bases], i, hb_C_d, s_C_d, hb_C_h,s_C_h, hb_C_a, s_C_a);
				fill_base(top, &base_bp[c_bases], i, hb_C_d, s_C_d, hb_C_h, s_C_h, hb_C_a, s_C_a);
				c_bases++;
				break;
			case 'g':
				//printf ("Found DG, %s\n", resname);
				allocate_base(top, &base_bp[c_bases], i, hb_G_d, s_G_d, hb_G_h,s_G_h, hb_G_a, s_G_a);
				fill_base(top, &base_bp[c_bases], i, hb_G_d, s_G_d, hb_G_h, s_G_h, hb_G_a, s_G_a);
				c_bases++;
				break;
			default:
				printf("Warning: could not identify %s as NA\n",res1);
				break;
		}
	}
	 

	for (i=0; i<c_bases; i++)
	{
		res1=*top.atoms.resinfo[base_bp[i].resnum].name;
		p1 = return_case(res1);
		for (j=i+1; j<(c_bases-1); j++)
		{
			res2=*top.atoms.resinfo[base_bp[j].resnum].name;
			p2 = return_case(res2);
			bp_type = return_bp_type(p1,p2);
			switch( bp_type )
			{
				case 1:
					// AT or TA
					hb_max[i][j] = 2;
					break;
				case 2:
					// GC or CG
					hb_max[i][j] = 3;
					break;
				case 0:
					hb_max[i][j] = 0;
					break;
				default: 
					hb_max[i][j] = 0;
					break;
			}
		}
	}
}

real compute_hb(int don, int h, int acc, rvec *x)
{
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

	rij[0] = x[h][0]-x[acc][0];
	rij[1] = x[h][1]-x[acc][1];
	rij[2] = x[h][2]-x[acc][2];

	mod_rij_d  = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
	rdist= mod_rij_d/0.32;
	if(rdist>0.99998 && rdist<1.00002) rdist = 0.99999;       // to keep the function continuos
	r6dist = pow(rdist, nn);
	num = 1.-r6dist;                                          // numerator
	r12dist = pow(rdist, mm);
	iden = 1./(1.-r12dist);
	// num/idem is the distance contribution to hbond
	d_contrib =  num*iden;

	rij[0] = x[don][0]-x[h][0];
	rij[1] = x[don][1]-x[h][1];
	rij[2] = x[don][2]-x[h][2];

	mod_rij  = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]) ;

	sij[0] = x[h][0]-x[acc][0];
	sij[1] = x[h][1]-x[acc][1];
	sij[2] = x[h][2]-x[acc][2];

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

	return(d_contrib*a_contrib);
}

real compute_all_hb(base *b, int num_base, rvec *x,real **hb_mat, int **hb_max, t_topology top)
{
	int i,j,k,l,m;
	real native_hb=0;
	int don=0, hyd=0, acc=0;
	real hbond;
	real total_bp=0;

	for(i=0; i< num_base; i++){
		for(j=i+1; j< (num_base-1); j++){
			hb_mat[i][j] = 0;
		}
	}
	for(i=0; i< num_base; i++){
		for(j=i+1; j< (num_base-1); j++){

			// first donors of i
			for(k=0; k<b[i].ndon; k++){
				don=b[i].don[k];
				hyd=b[i].h[k];
				for(l=0; l<b[j].nacc; l++){
					acc=b[j].acc[l];
					hbond=compute_hb(don,hyd,acc,x);
					if( hb_max[i][j] != 0 ){
						hb_mat[i][j] += (hbond/(real) hb_max[i][j]);
					//	hb_mat[i][j] += (hbond);
					}
				}
			}
			// then donors of j
			for(k=0; k<b[j].ndon; k++){
				don=b[j].don[k];
				hyd=b[j].h[k];
				for(l=0; l<b[i].nacc; l++){
					acc=b[i].acc[l];
					hbond=compute_hb(don,hyd,acc,x);
					if( hb_max[i][j] != 0 ){
						hb_mat[i][j] += (hbond/(real) hb_max[i][j]);
					//	hb_mat[i][j] += hbond;
					}
				}
			}

			total_bp += hb_mat[i][j];
		}
	}
	return total_bp;
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
	rvec       *xtop;
	matrix     box;
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
	int natoms;
	real t;
	t_trxstatus    *status;
	rvec *x;
	output_env_t    oenv;
	bIndex=TRUE;
	base *base_hb;
	int num_base=0;
	char *resname;
	gmx_rmpbc_t    gpbc = NULL;
	real **hb_mat;
	int **hb_max;

	//We should write out the density such that it can be re-analyzed afterwards...

	const char *desc[] = {
		"This g_program computes the number of native pairs hb. [BR]",
		"[BR]",
		"We use a smooth cut-off criteria based on the fulfillment of three criteria:[BR]",
		"[BR]",
		"------------------------------------------------------------------------------[BR]",
		" IMPORTANT!: So far PBC are not taken into account, so feed me well![BR]",
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
	};

#define NFILE asize(fnm)

	CopyRight(stderr,argv[0]);

	parse_common_args(&argc, argv, PCA_CAN_TIME | PCA_CAN_VIEW | PCA_BE_NICE,
			NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);

	int size_hb_T_d;
	int size_hb_T_h;
	int size_hb_T_a;
	int size_hb_A_d;
	int size_hb_A_h;
	int size_hb_A_a;
	int size_hb_C_d;
	int size_hb_C_h;
	int size_hb_C_a;
	int size_hb_G_d;
	int size_hb_G_h;
	int size_hb_G_a;
	int size_na;

	size_hb_T_d  = sizeof ( hb_T_d ) / sizeof (hb_T_d[0]) ;
	size_hb_T_h  = sizeof ( hb_T_h ) / sizeof (hb_T_h[0]) ;
	size_hb_T_a  = sizeof ( hb_T_a ) / sizeof (hb_T_a[0]) ;
	size_hb_A_d  = sizeof ( hb_A_d ) / sizeof (hb_A_d[0]) ;
	size_hb_A_h  = sizeof ( hb_A_h ) / sizeof (hb_A_h[0]) ;
	size_hb_A_a  = sizeof ( hb_A_a ) / sizeof (hb_A_a[0]) ;
	size_hb_C_d  = sizeof ( hb_C_d ) / sizeof (hb_C_d[0]) ;
	size_hb_C_h  = sizeof ( hb_C_h ) / sizeof (hb_C_h[0]) ;
	size_hb_C_a  = sizeof ( hb_C_a ) / sizeof (hb_C_a[0]) ;
	size_hb_G_d  = sizeof ( hb_G_d ) / sizeof (hb_G_d[0]) ;
	size_hb_G_h  = sizeof ( hb_G_h ) / sizeof (hb_G_h[0]) ;
	size_hb_G_a  = sizeof ( hb_G_a ) / sizeof (hb_G_a[0]) ;
	size_na  = sizeof ( na ) / sizeof (na[0]) ;


	// output
	out=ffopen(opt2fn("-o",NFILE,fnm),"w");
	//bIndex=ftp2bSet(efNDX,NFILE,fnm);
	bTPS = ftp2bSet(efTPS,NFILE,fnm) ;

	if ( bTPS ) {
		read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&xtop,NULL,box,TRUE);
		//sfree(xtop);
		top_atoms=&top.atoms;

	} else {
		printf("\nError: I need some sort of topology or pdb structure file\n\n");
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


	if (ePBC != epbcNONE) {
		snew(pbc,1);
	}else {
		pbc = NULL;	
	}

	natoms = read_first_x(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &t, &x, box);
	fflush(stdout);

    for (i=0; i<top.atoms.nres; i++)
	{
		resname = *(top.atoms.resinfo[i].name);
		if (in_list(na, size_na, resname) ){
			num_base++;
		}
	}
	snew(base_hb, num_base);
	snew(hb_mat, num_base);
	snew(hb_max, num_base);
	for (i=0; i<num_base; i++){
		snew(hb_mat[i],num_base);
		snew(hb_max[i],num_base);
	}
	//printf("\nAllocated %d bases\n", num_base);
	// Fill up the base_hb structure, determine max of hb per possible bp 
	populate_base(hb_T_d, size_hb_T_d, hb_T_h, size_hb_T_h, hb_T_a, size_hb_T_a, hb_A_d, 
			size_hb_A_d, hb_A_h, size_hb_A_d, hb_A_a, size_hb_A_a, hb_C_d, size_hb_C_d, 
			hb_C_h, size_hb_C_h, hb_C_a, size_hb_C_a, hb_G_d, size_hb_G_d, hb_G_h, 
			size_hb_G_h, hb_G_a, size_hb_G_a, natoms, top, base_hb, hb_max);
	printf("Done populating bases\n");

	if (pbc) {
		set_pbc(pbc,ePBC,box);
		gpbc = gmx_rmpbc_init(&top.idef, ePBC, natoms, box);
		gmx_rmpbc(gpbc, natoms, box, x);
	}

	real native_hb=0;
	do {
		native_hb=compute_all_hb(base_hb, num_base, x, hb_mat, hb_max, top, scheme);
		fprintf(out,"%f %f\n", t, native_hb);
	} while(read_next_x(oenv, status, &t, natoms, x, box));

	thanx(stderr);
	ffclose(out);
	return 0;

}

