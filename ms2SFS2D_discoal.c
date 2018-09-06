/******* ms2TwoSite2Popn.c ********


********************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "msGeneralStats.h"

#define LINEBUF 1000000



void usage();
void printSFS2D(int segsites, int nsam, int n1, int n2, char **list);
void tallySFS2D(int segsites, int nsam, int n1, int n2, char **list, int **jSFS);


int maxsites = 100000 ;
int sum=0;
int main(argc,argv)
	int argc;
char *argv[];
{
	int nsam, i,j,  howmany, **jSFS;
	char **list, **cmatrix(), line[LINEBUF+1]  ;
	FILE *fopen(), *pfin ;
	double *posit, theta ;
	int   segsites, count  , nadv, npops,n1, n2;
	 int nsites;
	char dum[20], astr[100] ;


/* read in first two lines of output  (parameters and seed) */
	pfin = stdin ;
	fgets( line, LINEBUF, pfin);
	sscanf(line," %s  %d %d %d -t %lf -p %d %d %d %s", dum,  &nsam, &howmany, &nsites, &theta, &npops,&n1,&n2, astr);
//	printf("dum: %s, theta: %lf npops: %d n1: %d rho: %f nsites: %d\n",dum,ss,npops,n1,rho,nsites);
	fgets( line, LINEBUF, pfin);
	
	if( argc > 1 ) { 
		nadv = atoi( argv[1] ) ; 
	}

	list = cmatrix(nsam,maxsites+1);
	jSFS = imatrix(n1+1,n2+1);
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
	//	printSFS2D(segsites,nsam,n1,n2,list);
		tallySFS2D(segsites,nsam,n1,n2,list,jSFS);
	}
	for(i=0;i<n1+1;i++){
		for(j=0;j<n2+1;j++){
			printf("%f\t",(float)jSFS[i][j]/sum);
		}
		printf("\n");
	}
	return(0);
}


void printSFS2D(int segsites, int nsam, int n1, int n2, char **list){
	int i, freq1, freq2;	 

	for(i=0; i < segsites; i++){
		freq1 = frequencySub('1', i, 0, n1,list);
		freq2 = frequencySub('1', i, n1, nsam,list);
		printf("%d\t%d\t%d\t%d\n",n1,n2,freq1,freq2);
	}
}

void tallySFS2D(int segsites, int nsam, int n1, int n2, char **list, int **jSFS){
	int i, freq1, freq2;	 

	for(i=0; i < segsites; i++){
		freq1 = frequencySub('1', i, 0, n1,list);
		freq2 = frequencySub('1', i, n1, nsam,list);
		jSFS[freq1][freq2]++;
		sum++;
	}
}


