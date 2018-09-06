/******* msMut.c ********
adds a backend finite site model to ms output
	********************************/

#include "msGeneralStats.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>


void printSFS(int segsites, int nsam, double *posit, int nsites, char **list);
void getParameters(int argc, char *argv[]);
void usage();

int maxsites = 1000 ;
int sites;


main(argc,argv)
	int argc;
char *argv[];
{
	int nsam, j ,nsites, i,  howmany  ;
	char **list, line[1001]  ;
	FILE *fopen(), *pfin ;
	double *posit,prob ;
	int   segsites, count,probflag  , nadv ;
	char dum[20], astr[100] ;
	int  segsub( int nsam, int segsites, char **list ) ;

//read in args
	getParameters(argc, argv);

/* read in first two lines of output  (parameters and seed) */
	pfin = stdin ;
	fgets( line, 1000, pfin);
	sscanf(line," %s  %d %d", dum,  &nsam, &howmany);
	fgets( line, 1000, pfin);

	if( argc > 1 ) { 
		nadv = atoi( argv[1] ) ; 
	}

	list = cmatrix(nsam,maxsites+1);
	posit = (double *)malloc( maxsites*sizeof( double ) ) ;

	count=0;
	probflag = 0 ;
	while( howmany-count++ ) {

/* read in a sample */
		do {
			fgets( line, 1000, pfin);
		}while ( line[0] != '/' );

		fscanf(pfin,"  segsites: %d", &segsites );
		if( segsites >= maxsites){
			maxsites = segsites + 10 ;
			posit = (double *)realloc( posit, maxsites*sizeof( double) ) ;
			biggerlist(nsam,maxsites, list) ;
		}
		if( segsites > 0) {
			fscanf(pfin," %s", astr);
			if( astr[1] == 'r' ){
				fscanf(pfin," %lf", &prob ) ;
				probflag = 1;
				fscanf(pfin," %*s");
			}
			for( i=0; i<segsites ; i++) fscanf(pfin," %lf",posit+i) ;
			for( i=0; i<nsam;i++) fscanf(pfin," %s", list[i] );
		}

/* analyse sample ( do stuff with segsites and list) */
		printf("//\n");
		printSFS(segsites,nsam,posit,sites,list);
	}
}



void getParameters(int argc, char *argv[]){

sites = atoi(argv[1]);	
}

void usage(){
  printf("usage:msMut nsites\n");
  exit(1);
}


void printSFS(int segsites, int nsam, double *posit, int nsites, char **list){
	int i, j, freq, nextSite, fdFlag;	 
	char **fsList;
	
	fsList = cmatrix(nsam,nsites);
	//initialize fsList
	for(i = 0; i<nsites;i++){
		for(j=0;j<nsam;j++){
			fsList[j][i]='0';
		}
	}
	
	//go though list, and translate to fsList
	for(i = 0; i<segsites;i++){
		nextSite = (int) floor(posit[i] * nsites);
		for(j=0;j<nsam;j++){
			fsList[j][nextSite]=list[j][i];
		}
	}	
	
	//print output	
	for(i = 0; i<nsites;i++){
		//check for fixedDiff
		if(fsList[0][i] == fsList[1][i]){
			fdFlag = 0;
		}
		else{
			fdFlag = 1;
		}
		
		freq = (((frequencySub('1', i, 1, nsam, fsList) > 0) && (frequencySub('0', i, 1, nsam, fsList) >0 )) ? 1:0);
		if(freq == 0 && fdFlag == 0){ printf("%d\t%d\t%d\n",i,0,nsam); }
		if(freq == 0 && fdFlag == 1){ printf("%d\t%d\t%d\n",i,1,nsam); }
		if(freq == 1 && fdFlag == 0){ printf("%d\t%d\t%d\n",i,2,nsam); }
		if(freq == 1 && fdFlag == 1){ printf("%d\t%d\t%d\n",i,3,nsam); }

	}
	
	for(j=0;j<nsam;j++){
		free(fsList[j]);
	}
	free(fsList);
}