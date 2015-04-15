#include <read_mult_t_inp.h>

void readfiles( const char *fn, f_info *fr, int *nfilesRet)
{
    char **filename = 0, tmp[255];
	char  line[255];
    int    nread, sizenow, i, block = 1;
    FILE  *fp;
	FILE  *f_melt;
	double xx,yy;
	int		j=0;

    fp      = ffopen(fn, "r");
    nread   = 0;
    sizenow = 0;
    while (fscanf(fp, "%s %lf %lf %lf %lf %lf %lf", tmp, 
				&fr[nread].conc, 
				&fr[nread].tm2, 
				&fr[nread].dh2, 
				&fr[nread].th2, 
				&fr[nread].bf, 
				&fr[nread].ef
			) != EOF)
	{
		fr[nread].tm2 += 273.15;  // we want it in Kelvins
		if (strlen(tmp) >= WHAM_MAXFILELEN)
		{
			gmx_fatal(FARGS, "Filename too long. Only %d characters allowed\n", WHAM_MAXFILELEN);
		}
		if (nread >= sizenow)
		{
			sizenow += block;
			srenew(fr[nread].fn, sizenow);
			for (i = sizenow-block; i < sizenow; i++)
			{
				snew(fr[i].fn, WHAM_MAXFILELEN);
			}
		}
		strcpy(fr[nread].fn, tmp);
		// now we read its data 
		if(tmp) {f_melt = ffopen(tmp,"r");}
		while(fgets(line, 255, f_melt) != NULL) {
			sscanf(line, "%lf %lf", &xx, &yy);
			if( xx > fr[nread].bf && xx < fr[nread].ef) {fr[nread].nts++; }
		}
		rewind(f_melt);
		snew(fr[nread].x[0],fr[nread].nts);
		snew(fr[nread].x[1],fr[nread].nts);
		snew(fr[nread].fx[0],fr[nread].nts);
		snew(fr[nread].fx[1],fr[nread].nts);
		//printf("Found %d frames in %s - l %f %f\n",fr[nread].nts, tmp, fr[nread].tm2, fr[nread].dh2);
		j=0;
		xx=0; yy=0;
		
		while(fgets(line, 255, f_melt) != NULL)
		{
			sscanf(line, "%lf %lf", &xx, &yy);
			if( xx > fr[nread].bf && xx < fr[nread].ef) {
				fr[nread].x[0][j]=xx; 
				fr[nread].x[1][j]=yy; 
				j++; 
			}
		}
		
		fclose(f_melt);

		nread++;
	}
   *nfilesRet    = nread;
   ffclose(fp);
}
