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
	int nsam, i,  howmany  ,nwins, siteIdx;
	char **list, **cmatrix(), line[LINEBUF+1]  ;
	FILE *fopen(), *pfin ;
	double *posit ;
	int   segsites, count  , nadv, iss,h ;
	double pi , th,  z, H, tajD,w, wins[50],max, min, temp_site;
	char dum[50], astr[100] ;



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
		nwins=9;
	printf("pi\tss\tthetaH\ttajD\tfayWuH\tHapCount\tZnS");
	/*printf("pi\tss\tthetaH\ttajD\tfayWuH\tHapCount");*/
	for(i=0;i<nwins;i++)printf("\tpiWin%d",i);
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


		iss = segSites(segsites,nsam,list);
		pi = nucdiv(nsam, segsites, list) ;
		th = thetah(nsam, segsites, list) ;
		h = nHaplotypes(segsites,nsam,list);
		H = th-pi;
		tajD = tajd(nsam,segsites,pi);
		z = ZnS( segsites,  nsam,  list);

		/*OmegaCenter*/
		/*Get snp that is nearest to our fixation*/
/*
			min = (abs(0.5 - posit[0]));
			for (i = 0; i < segsites; ++i){
				temp_site = posit[i];
					if (abs(0.5 - temp_site) <= min){
						siteIdx = i;
						min = abs(0.5-temp_site);
					}
			}
			w = omegaCenter(siteIdx, segsites, nsam, list);
*/		/*w = omegaMax(segsites, nsam,list);*/
		
		printf("%lf\t%d\t%lf\t%lf\t%lf\t%d\t%f", pi, iss, th ,tajD,H, h, z);
		//printf("%lf\t%d\t%lf\t%lf\t%lf\t%d", pi, iss, th ,tajD,H, h);
		
		//window stats
	
		nucdivWindow( nwins, posit, wins, nsam, segsites,list);
		//print normalized windows
		max=0.0;
		for( i=0; i<nwins ; i++){
			if (wins[i] > max){
				max = wins[i];
			}
		} 
		for( i=0; i<nwins ; i++){
			if (max == 0.0){ // Prevent "-nans" in output. 
 				printf("\t%f", 0.0);
			}
			else{
				printf("\t%f",wins[i]/max) ;
			} 
	}
	printf("\n");
 }
	return(0);
}
