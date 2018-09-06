/******* ms2SFSVector.c ********
converts ms output to SFS Vector
representation.

by convention fixed sites are 
treated as monomorphic
 
********************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "msGeneralStats.h"

#define LINEBUF 1000000



void usage();
int maxSampleSize = 5000;

int maxsites = 100000 ;

int main(argc,argv)
	int argc;
char *argv[];
{
	int nsam, i,  howmany  ;
	char **list, **cmatrix(), line[LINEBUF+1]  ;
	FILE *fopen(), *pfin ;
	double *posit, rho, mig ;
	int   segsites, count  , nadv, npops,n1, n2;
	 int nsites, derived_counts[maxSampleSize];
	double ss;
	char dum[20], astr[100] ;


/* read in first two lines of output  (parameters and seed) */
	pfin = stdin ;
	fgets( line, LINEBUF, pfin);
	sscanf(line," %s  %d %d  %s", dum,  &nsam, &howmany, astr);
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
	//	printf("n1 n2 p1 p2\n");
		getSiteFreqSpec(segsites,nsam,list,0,derived_counts);
		//print out the SFS header
		for(i=0;i<nsam;i++)printf("X_%d\t",i);
		printf("\n");
		for(i=0;i<nsam;i++)printf("%d\t",derived_counts[i]);
		printf("\n");
		printf("//\n");
	}
	return(0);
}




