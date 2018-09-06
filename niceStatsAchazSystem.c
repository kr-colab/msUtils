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
	int nsam, i,  howmany  ,nwins, siteIdx;
	char **list, **cmatrix(), line[LINEBUF+1]  ;
	FILE *fopen(), *pfin ;
	double *posit ;
	int   segsites, count  , nadv, iss, h, exponent1, exponent2;
	double pi , th,  z, H, tajD,w, wins[50],max, min, temp_site, thetaA, thetaHPi, tajDX, achazD, h1, h2, h12, thetaW;
        double winsH1[50], winsH2[50], winsH12[50];
	char dum[50], astr[100] ;
        double *harmonicSums;
        double ehh, rehh;
        int *haplotype_counts;


/* read in first two lines of output  (parameters and seed) */
	pfin = stdin ;
	fgets( line, LINEBUF, pfin);
	sscanf(line," %s  %d %d", dum,  &nsam, &howmany);
	fgets( line, LINEBUF, pfin);

	if( argc > 1 ) { 
		nadv = atoi( argv[1] ) ; 
	}

	list = cmatrix(nsam,maxsites+1);
	posit = (double *)malloc( maxsites*sizeof( double ) ) ;

	count=0;
	//print header
	nwins=5;


	printf("achazsD_1");
        for (exponent1 = 2; exponent1 < nsam; exponent1++){
                printf("\tachazsD_%d", exponent1);
        }

	printf("\n");
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
		}
		/* analyse sample ( do stuff with segsites and list) */

	        haplotype_counts = (int *)malloc( nsam*sizeof( int ) ) ;

		
		printf("%lf",achazThetaParabolicWeights(nsam,segsites,list,2,1.0));
		for (exponent1 = 2; exponent1 < nsam; exponent1++){
	                printf("\t%lf", achazThetaParabolicWeights(nsam,segsites,list,2,(double)exponent1));
	        }
		printf("\n");
		 
 }
	return(0);
}
