#include <stdlib.h>
#include <stdio.h>
#include "bitStuff.h"

int main(int argc,char *argv[]){
	int bit = 1000;
	
	printf("%d\n",bit);
	bit = setBit(bit,0);
	printf("%d\n",bit);
	bit = setBit(bit,1);
	printf("%d\n",bit);
	return(1);
}
