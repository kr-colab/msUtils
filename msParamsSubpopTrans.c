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

	xaRatio = 1.0;
	t1 = atof(argv[2]) * xaRatio;
	n1 = atof(argv[3]);
	t2 = atof(argv[4]) * xaRatio;
	thetaps = atof(argv[5]);
	rhops = atof(argv[6]);

	for(i = 0; i < n; i++){
		length = (data[i].chromEnd) - data[i].chromStart;
		theta = length *thetaps;
		rho = length*rhops;
		//adjust rho
		if (rho > 600){
			rho = 600;
		}

		//want theta, rho, length, trec, tbn, tAdmix, pAdmix,
		printf("%lf\t%lf\t%ld\t%0.12lf\t%0.12lf\t%0.12lf\t%0.12lf\t%0.12lf\n",theta,rho, length, t1,n1,t2, \
			thetaps, rhops);
	}
	
	
	return(0);
}

void usage(){
	printf("msParamsSubpopTrans bedFile <5 args>\n");
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
