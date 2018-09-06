/******* ms2TwoSite2Popn.c ********
converts ms output to list of all 
pairwise two-site comparisons for 2 popns.
The AFS in a two site 2 popn setting is a 6D
Matrix with the following entries
x = {p1, p2, x11, p3, p4, y11}
where p1 and p2 represent the frequency of
derived allele at locus one and locus two in popn 1,
X11 is the number of p1p2 haps in popn 1,
p3 and p4 are respective freqs in popn2,
and y11 is the number of p3,p4 haps in popn2

Seven column output: p1 p2 x11  p3 p4 y11 dist
********************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "msGeneralStats.h"

#define LINEBUF 1000000



void usage();
void printSFS2D(int segsites, int nsam, int n1, int n2, char **list, int *v);


int maxsites = 100000 ;
	 int nsites;
int main(argc,argv)
	int argc;
char *argv[];
{
	int nsam, i,  howmany  ;
	char **list, **cmatrix(), line[LINEBUF+1]  ;
	FILE *fopen(), *pfin ;
	double *posit, rho, mig ;
	int   segsites, count  , nadv, npops,n1, n2, *positSites;

	double ss, t;
	char dum[20], astr[100] ;


/* read in first two lines of output  (parameters and seed) */
	pfin = stdin ;
	fgets( line, LINEBUF, pfin);
	sscanf(line," %s  %d %d %d -t %lf -r %lf -p %d %d %d %s", dum,  &nsam, &howmany, &nsites, &t, &rho ,&npops,&n1,&n2, astr);
//	printf("dum: %s, theta: %lf npops: %d n1: %d rho: %f nsites: %d\n",dum,ss,npops,n1,rho,nsites);
	fgets( line, LINEBUF, pfin);
	
	if( argc > 1 ) { 
		nadv = atoi( argv[1] ) ; 
	}

	list = cmatrix(nsam,maxsites+1);
	posit = (double *)malloc( maxsites*sizeof( double ) ) ;
	positSites = (int *)malloc(nsites*sizeof(int));
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
			for( i=0; i<segsites ; i++) positSites[i] = floor(posit[i]*nsites);
			for( i=0; i<nsam;i++) fscanf(pfin," %s", list[i] );
		}
		/* analyse sample ( do stuff with segsites and list) */
	//	printf("n1 n2 p1 p2\n");
		printSFS2D(segsites,nsam,n1,n2,list,positSites);
		printf("//\n");
	}
	return(0);
}


void printSFS2D(int segsites, int nsam, int n1, int n2, char **list, int *posit){
	int i,j, freq1, freq2;	 
	char allele;

	j=0;
	for(i=0; i < segsites; i++){
		while(j<posit[i]){
			printf("%d\t%d\t%d\t%d\t%d\n",j,n1,n2,0,0);
			j++;
		}
		j++;
		freq1 = frequencySub('1', i, 0, n1,list);
		freq2 = frequencySub('1', i, n1, nsam,list);
		printf("%d\t%d\t%d\t%d\t%d\n",posit[i],n1,n2,freq1,freq2);
	}
	while(j<nsites){
		printf("%d\t%d\t%d\t%d\t%d\n",j,n1,n2,0,0);
		j++;
	}
}



