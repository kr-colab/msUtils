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
int maxsites = 100000 ;
void usage();

int main(argc,argv)
	int argc;
char *argv[];
{
	int nsam, i,  howmany  ;
	char **list, **cmatrix(), line[LINEBUF+1]  ;
	FILE *fopen(), *pfin ;
	double *posit ;
	int   segsites, count  , n1, n2, iss,h ;
	double pi , th,  z, f, snn,dxy, dxy_min, dxy_mean, H, tajD;
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
	else{
		usage();
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


		iss = segSitesSub(segsites,nsam,0,n1,list);
		pi = nucdivSub(nsam,segsites,0,n1,list);
		th = thetahSub(nsam, segsites,0,n1, list) ;
		h = nHaplotypesSub(segsites,nsam,0,n1,list);
		z = ZnSSub( segsites,  nsam, 0,n1, list);
		f = fst2Subs(segsites,nsam,0,n1,n1,n1+n2,list);
		H = th-pi;
		tajD = tajd(nsam,iss,pi);
		//printf("%f\t%f\t%f\n",nucdiv(nsam,segsites,list),nucdivSub(nsam,segsites,0,nsam,list),);
		printf("pi1:\t%lf\tss1:\t%d\tthetaH1:\t%lf\ttajd1:\t%f\tH1:\t%f\tHapCount1:\t%d\tZnS1:\t%f\t",pi, iss, th,tajD,H , h, z) ;
		iss = segSitesSub(segsites,nsam,n1,nsam,list);
		pi = nucdivSub(nsam,segsites,n1,nsam,list);
		th = thetahSub(nsam, segsites,n1,nsam, list) ;
		h = nHaplotypesSub(segsites,nsam,n1,nsam,list);
		z = ZnSSub( segsites,  nsam, n1,nsam, list);
		snn = Snn(segsites,nsam,n1,n2,list);
		dxy=Dxy(segsites,nsam,n1,n2,list);
		dxy_min = Dxy_min(segsites,nsam,n1,n2,list);
		dxy_mean = Dxy_mean(segsites,nsam,n1,n2,list);
		H = th-pi;
		tajD = tajd(nsam,iss,pi);
		printf("pi2:\t%lf\tss2:\t%d\tthetaH2:\t%lf\ttajd2:\t%f\tH2:\t%f\tHapCount2:\t%d\tZnS2:\t%f\tFst:\t%f\tsnn:\t%f\tdxy:\t%f\tdxy_mean:\t%f\tdxy_min:\t%f\n", pi, iss, th,tajD,H , h, z,f,snn,dxy,dxy_mean,dxy_min) ;


	}
	return(0);
}

void usage(){
	printf("maskedStatsSubpop n1 n2\n");
	printf("returns analysis of Hudson style output assuming two subpops of size n1 and n2\n");
	exit(1);
}

