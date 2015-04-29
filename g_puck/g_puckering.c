// GPC 29/04/2015 g_puckering v 1.0
// 
// Computes puckering angle and amplitude for ribose DNA/RNA rings
// Computes average of the puckering angle
// Reads the residues from a list in a file
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
#include <gromacs/bondf.h>
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

t_atoms  *top_atoms=NULL;

char *nucleic[]= {  "ABU", "ACE", "ACE", "AIB", "ALA", "ARG", "ARGN", "ASH", "ASN", "ASN1", "ASP", "ASP1", "ASPH", "CALA", "CARG", "CASN", "CASP", "CCYN", "CCYX", "CGLN", "CGLU", "CGLY", "CHID", "CHIE", "CHIP", "CILE", "CLEU", "CLYP", "CMET", "CPHE", "CPRO", "CSER", "CTHR", "CTRP", "CTYR", "CVAL", "CYM", "CYN", "CYS", "CYS1", "CYS2", "CYSH", "CYX", "DA", "DA3", "DA5", "DAB", "DALA", "DAN", "DC", "DC3", "DC5", "DCN", "DG", "DG3", "DG5", "DGN", "DT", "DT3", "DT5", "DTN", "GLH", "GLN", "GLU", "GLUH", "GLY", "HID", "HIE", "HIP", "HIS", "HIS1", "HISA", "HISB", "HISH", "HYP", "ILE", "LEU", "LYN", "LYP", "LYS", "LYSH", "MELEU", "MET", "MEVAL", "NAC", "NALA", "NARG", "NASN", "NASP", "NCYN", "NCYX", "NGLN", "NGLU", "NGLY", "NH2", "NHE", "NHID", "NHIE", "NHIP", "NILE", "NLEU", "NLYP", "NME", "NMET", "NPHE", "NPRO", "NSER", "NTHR", "NTRP", "NTYR", "NVAL", "ORN", "PHE", "PHEH", "PHEU", "PHL", "PRO", "RA", "RA3", "RA5", "RAN", "RC", "RC3", "RC5", "RCN", "RG", "RG3", "RG5", "RGN", "RU", "RU3", "RU5", "RUN", "SER", "THR", "TRP", "TRPH", "TRPU", "TYR", "TYRH", "TYRU", "VAL", "NA5", "NA3", "NA", "NG5", "NG3", "NG", "NC5", "NC3", "NC", "NU5", "NU3", "NU", "NT5", "NT3", "NT", "FA5", "FA3", "FA", "FG5", "FG3", "FG", "FC5", "FC3", "FC", "FU5", "FU3", "FU", "FT5", "FT3", "FT", "A","T","U","C","G", "AMO", "GMO", "TMO","CMO","PMO" "DP5", "DP", "DP3" };

char *pur[]={ "N1","C2","N3","C4","C5","C6","N7","C8","N9"};
char *pyr[]={ "N1","C2","N3","C4","C5","C6"};

char *bpur[]={"DG5", "DG", "DG3","DA", "DA3", "DA5", "RG5", "RG", "RG3","RA", "RA3", "RA5", "A","G", "AMO", "GMO", "A5E", "A3E"  };
char *bpyr[]={"DC", "DC3", "DC5", "DT", "DT3", "DT5","RC", "RC3", "RC5", "RU", "RU3", "RU5", "C","T", "U", "TMO", "CMO", "PMO", "DP", "DP5", "DP3", "T5E" };

typedef struct{
	int  res;
	int  o4;
	int  c1;
	int  c2;
	int  c3;
	int  c4;
	real  p;
	real  sp;
	real  cp;
	real  av_p;
	real  amp;
} sugar;

void seach_rings(int natoms,t_topology top, int nres, sugar *ring)
{
	int i,j;
	int sec;
	char *name;
	sec=0; 
	for (i=0; i<nres; i++){
		for (j=0; j<natoms; j++){
			if(top.atoms.resinfo[top.atoms.atom[j].resind].nr == ring[i].res) {
				name = *top.atoms.atomname[j] ; 
				if (strcmp("O4'",name) == 0) { 
					ring[i].o4 = j ;
					sec++;
				} else if (strcmp("C1'",name) == 0) {
					ring[i].c1 = j ;
					sec++;
				} else if (strcmp("C2'",name) == 0) {
					ring[i].c2 = j ;
					sec++;
				} else if (strcmp("C3'",name) == 0) {
					ring[i].c3 = j ;
					sec++;
				} else if (strcmp("C4'",name) == 0) {
					ring[i].c4 = j ;
					sec++;
				}
			}
		}
		if (sec!=5) { gmx_fatal(FARGS, "Could not find all 5 atoms for the ring in %d\n",ring[i].res); sec=0;}
		sec=0;
	}

}


static inline real dihedral(const rvec xi,const rvec xj,const rvec xk,const rvec xl)

{
  real ipr,phi,cos_phi,sign;
  rvec r_ij,r_kj,r_kl,m,n;

  rvec_sub(xi,xj,r_ij);
  rvec_sub(xj,xk,r_kj);
  rvec_sub(xk,xl,r_kl);

  cprod(r_ij,r_kj,m);
  cprod(r_kj,r_kl,n);
  cos_phi=cos_angle(m,n);
  if(cos_phi < -1.) cos_phi = -1.;
  else if(cos_phi > 1.) cos_phi = 1.;
  phi=acos(cos_phi);
  ipr=iprod(r_ij,n);
  sign=(ipr>0.0)?-1.0:1.0;
  phi=sign*phi;

  return phi;
}
void compute_p (int natoms,t_topology top,int nres,sugar *ring,int  nframes,FILE *f_back, FILE *f_amp, rvec *x, t_pbc *pbc)
{

	int i,j,k;
	rvec *at;
	real *v;
	real A, B, tm, c, s, t;
    int  t1, t2, t3;
	bool	bCurves;
	bCurves = FALSE;
    rvec r_ij, r_kj, r_kl, m, n;
	real sign;
	snew(v,5);
	snew(at,5);
	fprintf(f_back,"%d ", nframes);
	fprintf(f_amp,"%d ", nframes);
	for (i=0; i<nres; i++){
		copy_rvec(x[ring[i].o4],at[0]);	
		copy_rvec(x[ring[i].c1],at[1]);	
		copy_rvec(x[ring[i].c2],at[2]);	
		copy_rvec(x[ring[i].c3],at[3]);	
		copy_rvec(x[ring[i].c4],at[4]);	
		//v0 c4'-o4'-c1'-c2'
		//v1 o4'-c1'-c2'-c3'
		//v2 c1'-c2'-c3'-c4'
		//v3 c2'-c3'-c4'-o4'
		//v4 c3'-c4'-o4'-c1'
		//v[0] = dihedral(at[4],at[0],at[1],at[2]);
		v[0] = dih_angle(at[4],at[0],at[1],at[2],pbc, r_ij, r_kj, r_kl, m, n, &sign, &t1, &t2, &t3);
		v[1] = dih_angle(at[0],at[1],at[2],at[3],pbc, r_ij, r_kj, r_kl, m, n, &sign, &t1, &t2, &t3);
		v[2] = dih_angle(at[1],at[2],at[3],at[4],pbc, r_ij, r_kj, r_kl, m, n, &sign, &t1, &t2, &t3);
		v[3] = dih_angle(at[2],at[3],at[4],at[0],pbc, r_ij, r_kj, r_kl, m, n, &sign, &t1, &t2, &t3);
		v[4] = dih_angle(at[3],at[4],at[5],at[1],pbc, r_ij, r_kj, r_kl, m, n, &sign, &t1, &t2, &t3);
		//printf("V0 is %f V1 %f V2 %f V3 %f V4 %f\n", v[0],v[1],v[2],v[3],v[4]);
		if(bCurves){
			A = 0;
			B = 0;
			for (j=0; j<5; j++){
				t = 0.8 * PI * (j-1);
				A += v[j] * cos(t);
				B += v[j] * sin(t);
			}
			A *= 0.4;
			B *= -0.4;
			tm = sqrt(A*A + B*B);
			ring[i].amp = tm;
			c = A/tm ; 
			s = B/tm ; 
			ring[i].p = atan2(s,c);
		}else{
			A = sin(PI/5) + sin(PI/2.5);
			//printf("V0 is %f V1 %f V2 %f V3 %f V4 %f\n", v[0],v[1],v[2],v[3],v[4]);
			ring[i].p = atan2(v[4] + v[1] - v[3] - v[0], 2.0 * v[2] * A); 
			ring[i].amp = v[2] / cos(ring[i].p); 
		}

		ring[i].sp += sin(ring[i].p);
		ring[i].cp += cos(ring[i].p);
		fprintf(f_back,"%f ", ring[i].p);
		fprintf(f_amp,"%f ", ring[i].amp);
	}

	fprintf(f_back,"\n");
	fprintf(f_amp,"\n");
	sfree(at);
	sfree(v);
}

bool is_in_my_list(char *list[], int size, char *resnm )
{
  /* gives true if resnm occurs in the list of names provided */
  int i;

  for(i=0; (i<size); i++)
    //if (strcasecmp(nucleic[i],resnm) == 0)
    if (strcmp(list[i],resnm) == 0)
      return TRUE;

  return FALSE;
}



int main(int argc, char *argv[])
{
	const char         *desc[] = {
		"Read a pdb + trajectory, compute pseudorotation angle[BR]"
		"for the ribose rings of the residues listed in the index file.[BR]"
		"Output the time dep. pseudorotation angle and amplitudes, and the residue[BR]"
		"averaged pseudorotation angle in differents files."
	};

	FILE				*f_res;
	FILE				*f_back;
	FILE				*f_amp;
	FILE				*f_prof;
	static bool	bVerbose=FALSE;

	/* Command-line arguments */
	t_pargs          pa[] = {
		{"-v", TRUE, etBOOL, {&bVerbose},
			"Be verbose, output correlations for each replica to std output."},
	};

	/* Output files */
	t_filenm            fnm[] = {
		{ efTRX, "-f", NULL, ffREAD },
		{ efTPS, NULL, NULL, ffREAD },
		{ efDAT, "-r", "res_pucker", ffREAD },
		{ efXVG, "-o", "puckering", ffWRITE  },
		{ efXVG, "-pr", "puck_profile", ffWRITE  },
		{ efXVG, "-a", "amplitude", ffWRITE  },
	};

#define NFILE asize(fnm)
	output_env_t	oenv;
	t_tpxheader		header;
	matrix			box;
	t_topology      top;
	t_trxstatus		*status;
	t_fileio		*trxout;
	real            t;
	rvec			*x,*xtop;
	int				natoms,i,j,kkk;
	int				ePBC;
	int				nframes=0;
	int				nres;
	int				n_aa;
	int				step=1;
	int				nfr;
	char			buf[256];
	const char		*out_file=NULL; 
	const char		*resfile=NULL; 
	const char		*puckfile=NULL; 
	const char		*ampfile=NULL; 
	const char		*profile=NULL; 
	char			line[255];
	char       title[STRLEN];
	t_pbc 			*pbc;
	snew(pbc, 1);

	int npargs;
	npargs=asize(pa);
	CopyRight(stderr, argv[0]);
	parse_common_args(&argc, argv, PCA_CAN_TIME | PCA_BE_NICE,
			NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);

    if ( ftp2bSet(efTPS,NFILE,fnm) ) {
        read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&xtop,NULL,box,TRUE);
        //sfree(xtop);
        top_atoms=&top.atoms;

    } else {
        printf("\nSTOP: I need some sort of topology or pdb structure file\n\n");
        exit(1);
    }

    if (pbc)
	{
		set_pbc(pbc, -1, box);
	}

	// Read RDC file
	resfile = opt2fn("-r",NFILE,fnm);
	if(resfile) {f_res = ffopen(resfile,"r");}
	nres=0;
	while(fgets(line, 255, f_res) != NULL)
	{
		nres++;
	}

	if (bVerbose) {
		printf("There are %d residues from %s\n",nres,resfile);
	}
	rewind(f_res);

	if ( nres < 1  ) {gmx_warning("The number of residues should at least be one");}

	sugar *ring;
	snew(ring,nres);
	for (i=0; i< nres ; i++){
		ring[i].p = 0;
		ring[i].sp = 0;
		ring[i].cp = 0;
		ring[i].av_p = 0;
	}

	// Crude parsing, no comments or extra lines allowed
	j=0;
	while(fgets(line, 255, f_res) != NULL)
	{
		sscanf(line, "%d ", &ring[j].res);
		j++;
	}

	int  size_of_aa_record;
	size_of_aa_record= sizeof ( nucleic ) / sizeof (nucleic[0]) ;
	for (i=0; i< nres; i++){
		if(is_in_my_list( nucleic, size_of_aa_record, *top.atoms.resinfo[ring[i].res -1 ].name) == FALSE){
			gmx_fatal(FARGS,"Residue %d is not a nucleic acid (%s)", ring[i].res,  *top.atoms.resinfo[ring[i].res -1 ].name);
		}
	}

	// Read atoms from pdb file
	natoms = read_first_x(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &t, &x, box);
	// find rings
	seach_rings(natoms, top, nres, ring);

	// Output computed vs experimental RDCs
	puckfile = opt2fn("-o",NFILE,fnm);
	ampfile = opt2fn("-a",NFILE,fnm);
	profile = opt2fn("-pr",NFILE,fnm);
	if(puckfile && ampfile && profile) {
		f_back = ffopen(puckfile,"w");
		f_amp = ffopen(ampfile,"w");
		f_prof = ffopen(profile,"w");
	}

	nfr=0;
	do {

		compute_p(natoms, top, nres, ring, nfr, f_back, f_amp,  x, pbc);
		nfr++;

	} while (read_next_x(oenv, status, &t, natoms, x, box));
	printf("Read %d frames\n", nframes);

	if(profile) {
		for (i=0; i< nres ; i++){
			ring[i].av_p = atan2(ring[i].sp/(real)nfr, ring[i].cp/(real)nfr);
			if ( ring[i].av_p < 0 ) { ring[i].av_p = 2*PI + ring[i].av_p ;}
			fprintf(f_prof,"%d %f\n", ring[i].res, ring[i].av_p*180.0/PI);
		}
	}

	if(puckfile && ampfile && profile) {
		fclose(f_back);
		fclose(f_amp);
		fclose(f_prof);
	}


	return 0;

}

