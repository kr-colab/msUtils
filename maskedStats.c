/******* maskedStats.c ********
for calculating sample stats from MS output 
after it has been filtered by msMask
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
	int nsam, i,  howmany  ;
	char **list, **cmatrix(), line[LINEBUF+1]  ;
	FILE *fopen(), *pfin ;
	double *posit ;
	int   segsites, count  , nadv, iss,h ;
	double pi , th,  z;
	char dum[20], astr[100] ;



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


		iss = segSites(segsites,nsam,list);
		pi = nucdiv(nsam, segsites, list) ;
		th = thetah(nsam, segsites, list) ;
		h = nHaplotypes(segsites,nsam,list);
		z = ZnS( segsites,  nsam,  list);
		//printf("%f\t%f\t%f\n",nucdiv(nsam,segsites,list),nucdivSub(nsam,segsites,0,nsam,list),fst2Subs(segsites,nsam,0,5,5,10,list));
		printf("pi:\t%lf\tss:\t%d\tthetaH:\t%lf\tHapCount:\t%d\tZnS:\t%f\n", pi, iss, th , h, z) ;


	}
	return(0);
}




