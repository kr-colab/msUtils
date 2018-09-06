/******* slideFST.c ********

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
	int nsam, i,  howmany, siteIdx;
	char **list, **cmatrix(), line[LINEBUF+1]  ;
	FILE *fopen(), *pfin ;
	double *posit ;
	int 	nwins=200;
	int   segsites, count  ,n1,n2 ;
	double  wins[nwins],max, min, temp_site;
	char dum[20], astr[100] ;



/* read in first two lines of output  (parameters and seed) */
	pfin = stdin ;
	fgets( line, LINEBUF, pfin);
	sscanf(line," %s  %d %d", dum,  &nsam, &howmany);
	fgets( line, LINEBUF, pfin);

	if( argc > 1 ) { 
			n1 = atoi( argv[1] ) ; 
			n2 = atoi( argv[2] ) ;
	}

	list = cmatrix(nsam,maxsites+1);
	posit = (double *)malloc( maxsites*sizeof( double ) ) ;

	count=0;
	//print header

	printf("fstWin0");
	for(i=1;i<nwins;i++)printf("\tpiWin%d",i);
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
		
		//window stats		
		fst2SubsWindow( nwins, posit, wins, segsites, nsam, 0,n1,n1,n1+n2,list);
		//print  windows
		for( i=0; i<nwins ; i++) printf("%f\t",wins[i]) ;
		printf("\n");

	}
	return(0);
}




