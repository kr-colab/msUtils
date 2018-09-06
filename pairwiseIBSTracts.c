/******* niceStats.c ********
for calculating sample stats from MS output 
********************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "msGeneralStats.h"

#define LINEBUF 1000000



void usage();



int maxsites = 100000 ;

int main(argc,argv)
	int argc;
char *argv[];
{
	int nsam, i, j, k, howmany;
	char **list, **cmatrix(), line[LINEBUF+1]  ;
	FILE *fopen(), *pfin ;
	double *posit ;
	int   segsites, count;
	double tmpLen;

	char dum[50], astr[100] ;
        double start;
	double sites;
	int bins = 10;
	double binWidth = 1.0 / ((float) bins);
	int hist[bins+1];
	float sumComp = 0;
	
	
	if(argc < 1){
		usage();
	}
/* read in first two lines of output  (parameters and seed) */
	pfin = stdin ;
	fgets( line, LINEBUF, pfin);
	sscanf(line," %s  %d %d", dum,  &nsam, &howmany);
	fgets( line, LINEBUF, pfin);


	sites = atof(argv[1]);	
	list = cmatrix(nsam,maxsites+1);
	posit = (double *)malloc( maxsites*sizeof( double ) ) ;

	count=0;
	

	while( howmany-count++ ) {

/* read in a sample */
		do {
			fgets( line, LINEBUF, pfin);
		}while ( line[0] != '/' );

		fscanf(pfin,"  segsites: %d", &segsites );
		if( segsites >= maxsites){
			maxsites = segsites + 10 ;
			posit = (double *)realloc( posit, maxsites*sizeof( double) ) ;
			biggerlist(nsam,maxsites, list) ;
		}
		if( segsites > 0) {
			fscanf(pfin," %s", astr);

			for( i=0; i<segsites ; i++) fscanf(pfin," %lf",posit+i) ;
			for( i=0; i<nsam;i++) fscanf(pfin," %s", list[i] );
			for( i=0; i<bins+1;i++) hist[i]=0;
			sumComp=0;
		}
		/* analyse sample ( do stuff with segsites and list) */
		//first iterate through pairwise comps
		for( i=0; i<nsam-1;i++){
			for(j=i+1;j<nsam;j++){
				start=0.0;
				//now iterate across sites
				for(k=0;k<segsites;k++){
					if(list[i][k] != list[j][k]){
						//tmpLen = ((int) round((posit[k]-start) * sites))-1;
						tmpLen = posit[k]-start;
						//printf("%f\n",round(tmpLen/binWidth));
						hist[(int)round(tmpLen/binWidth)]+= 1;
						sumComp +=1;
						//if(tmpLen > 0)printf("%d\n",tmpLen);
						start = posit[k];
					}
				}
				//tmpLen =(int) round((1.0-start) * sites);
				tmpLen = 1.0-start;
				//printf("%f\n",round(tmpLen/binWidth));
				hist[(int)round(tmpLen/binWidth)]+= 1;
				sumComp +=1;
				//if(tmpLen > 0)printf("%d\n",tmpLen);
			}
		}

		for (i=0;i<bins+1;i++) printf("%lf\t",hist[i]/sumComp);
		printf("\n"); 
 
 	}
	return(0);
}

void usage(){
	printf("pairwiseIBSTracts totalSites -- returns distrb. of all pairwise IBS tract lengths\n");
	exit(6);
}