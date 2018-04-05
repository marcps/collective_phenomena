#include <stdio.h>
#include <stdlib.h>

//calls the Fortran SUBROUTINE RCARIN(IJKL,RVEC,LENV)
void rcarin_(int *,float *, int *);
//SUBROUTINE RCARRY(RVEC,LENV)
void rcarry_(float *, int *);

int main(void){
	int nrand=1000;
	float rrand[nrand+24];
	int i,seed;
	
	//Initialization
	for(i=0;i<nrand+24;i++){
		rrand[i]=0.;
	}
	rcarin_(&seed,rrand,&nrand);
	rcarry_(rrand,&nrand);
	
	for(i=0;i<nrand;i++){
		printf("%f\n",rrand[i]);
	}
	return 0;
}
