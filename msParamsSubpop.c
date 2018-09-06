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
	double t1, t2, n1, n2, theta, rho, rhops, thetaps,nAfr;
	double pAdmix, tAdmix,xaRatio;
	
	
	if(argc < 2){
		usage();
		exit(1);
	}
	seed1 = (long) abs(devrand() % 2147483399);
	seed2 = (long) abs(devrand() % 2147483399);
	setall(seed1, seed2);
	n = bedFileImport3(argv[1],data);
	if(argc < 3){
		t1 = genunf(0.0,0.2);
		t2 = genunf(t1,0.8);
		n1 = genunf(0.0,0.2);
	//	n2 = genunf(n1,1.0);
		tAdmix = genunf(0,t2);
		pAdmix = genunf(0,0.75);
//		nAfr = genunf(0,1.5);
		thetaps = genunf(0.00000001,0.03);
		rhops = genunf(0,0.3);
	}
	else{
		xaRatio = atof(argv[7]);
		tAdmix = atof(argv[2]) * xaRatio;
		pAdmix = atof(argv[3]);
		t1 = atof(argv[4]) * xaRatio;
		n1 = atof(argv[5]);
		t2 = atof(argv[6]) * xaRatio;
		thetaps = genunf(0.00000001,0.03);
		rhops = genunf(0,0.3);
	}
	for(i = 0; i < n; i++){
		length = (data[i].chromEnd) - data[i].chromStart;
		theta = length *thetaps;
		rho = length*rhops;
		//adjust rho
		if (rho > 600){
			rho = 600;
		}

		//want theta, rho, length, trec, tbn, tAdmix, pAdmix,
		printf("%lf\t%lf\t%ld\t%0.12lf\t%0.12lf\t%0.12lf\t%0.12lf\t%0.12lf\t%0.12lf\t%0.12lf\n",theta,rho, length, tAdmix, pAdmix, t1,n1,t2, \
			thetaps, rhops);
	}
	
	
	return(0);
}


void usage(){
	printf("msParamsSubpop bedFile <optional 7 params> <x/auto ratio>\n");
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
