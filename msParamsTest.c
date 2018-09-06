/*msParamsSubpop.c **************
/ makes a params file for 2
/ popn scenario
***************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fcntl.h>
#include <unistd.h>
#include "../coalLib/ranlib.h"
#include "../pgLib/bedFile.h"

void usage();
unsigned int devrand(void);

int main(int argc, char *argv[]){
	int i,n;
	long int length, seed1, seed2;
	struct bedEl data[50000];
	double t1, t2, n1, n2, theta, rho, rhoMax, thetaMax,nAfr;
	double pAdmix, tAdmix;
	
	
	if(argc < 2){
		usage();
		exit(1);
	}
	seed1 = (long) abs(devrand() % 2147483399);
	seed2 = (long) abs(devrand() % 2147483399);
	setall(seed1, seed2);
	
	n = bedFileImport3(argv[1],data);
	t1 = genunf(0.0,0.2);
	t2 = genunf(0.0,0.8);
	n1 = genunf(0.0,0.2);
	n2 = genunf(n1,1.0);
	tAdmix = genunf(0,t2);
	pAdmix = genunf(0,0.75);
	nAfr = genunf(0.8,3.0);
	for(i = 0; i < n; i++){
		length = (data[i].chromEnd) - data[i].chromStart;
		thetaMax = length * 0.02 * 10.0;
		theta = genunf(0.001,thetaMax);
		rhoMax = thetaMax * 5.0;
		rho = genunf(0,rhoMax);
		//adjust rho
		if (rho > 600){
			rho = 600;
		}

		//want theta, rho, length, trec, tbn, tAdmix, pAdmix,
		printf("%lf\t%lf\t%ld\t%0.12lf\t%0.12lf\n",theta,rho, length, nAfr,t2);
	}
	
	
	return(0);
}


void usage(){
	printf("msParams bedFile\n");
}

/* used for getting random number seeds */
unsigned int devrand(void) {
	int fn; 
	unsigned int r; 

	fn = open("/dev/urandom", O_RDONLY); 
	if (fn == -1) 
		exit(-1); /* Failed! */ 
	if (read(fn, &r, 4) != 4) 
		exit(-1); /* Failed! */ 
	close(fn); 
	return r;
}
