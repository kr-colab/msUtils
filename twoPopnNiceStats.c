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
	int   segsites, count  , n1, n2, iss,h,nnm1 ;
	double pi , th,  z, f, snn,dxy, dxy_min, dxy_mean, H, tajD, *dxyVec;
	char dum[20], astr[100] ;



/* read in first two lines of output  (parameters and seed) */
	pfin = stdin ;
	fgets( line, LINEBUF, pfin);
	sscanf(line," %s  %d %d", dum,  &nsam, &howmany);
	fgets( line, LINEBUF, pfin);
	
	nnm1=0;
	if( argc > 1 ) { 
		n1 = atoi( argv[1] ) ; 
		n2 = atoi( argv[2] ) ; 
		nnm1 = n1 * n2;
	}
	else{
		usage();
	}

	list = cmatrix(nsam,maxsites+1);
	posit = (double *)malloc( maxsites*sizeof( double ) ) ;
	dxyVec = (double *)malloc( nnm1*sizeof( double ) ) ;
	count=0;
	//print header
//	printf("pi1\tss1\tthetaH1\ttajd1\tH1\tHapCount1\tZnS1\t") ;
//	printf("pi2\tss2\tthetaH2\ttajd2\tH2\tHapCount2\tZnS2\t") ;
//	printf("Fst\tsnn\tdxy\tdxy_mean\tdxy_min\n");
	
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

		
		//now get polymorphism numbers for each and print
		iss = segSitesSub(segsites,nsam,0,n1,list);
		pi = nucdivSub(nsam,segsites,0,n1,list);
		th = thetahSub(nsam, segsites,0,n1, list) ;
		H = th-pi;
		h = nHaplotypesSub(segsites,nsam,0,n1,list);
		z = ZnSSub( segsites,  nsam, 0,n1, list);
		tajD = tajd(n1,iss,pi);
		printf("%lf\t%d\t%lf\t%lf\t%lf\t%d\t%lf\t",pi, iss, th,tajD,H , h, z) ;
		iss = segSitesSub(segsites,nsam,n1,nsam,list);
		pi = nucdivSub(nsam,segsites,n1,nsam,list);
		th = thetahSub(nsam, segsites,n1,nsam, list) ;
		H = th-pi;
		h = nHaplotypesSub(segsites,nsam,n1,nsam,list);
		z = ZnSSub( segsites,  nsam, n1,nsam, list);
		tajD = tajd(n2,iss,pi);
		printf("%lf\t%d\t%lf\t%lf\t%lf\t%d\t%lf\t",pi, iss, th,tajD,H , h, z) ;
		
		//subdivision stats
		f = fst2Subs(segsites,nsam,0,n1,n1,n1+n2,list);
		snn = Snn(segsites,nsam,n1,n2,list);
		printf("%lf\t%lf\t",f,snn) ;
		
		dxy=Dxy(segsites,nsam,n1,n2,list);
		dxy_min = Dxy_min(segsites,nsam,n1,n2,list);
		dxy_mean = Dxy_mean(segsites,nsam,n1,n2,list);
		//printf("%lf\t",dxy_min/dxy_mean);
		Dxy_vector(segsites,nsam, n1, n2, dxyVec,list);
		for(i=0;i<nnm1-1;i++)printf("%lf\t",dxyVec[i]/dxy_mean);
		printf("%lf\n",dxyVec[i]/dxy_mean);
		//printf("\n");

	}
	return(0);
}

void usage(){
	printf("twoPopnNiceStats n1 n2\n");
	printf("returns analysis of Hudson style output assuming two subpops of size n1 and n2\n");
	exit(1);
}

