#include <stdlib.h>
#include <stdio.h>


//calls the Fortran SUBROUTINE RCARIN(IJKL,RVEC,LENV)
//rcarin_(&seed,rrand,&nrand);
//rcarry_(rrand,&nrand);
void rcarin_(int *,float *, int *);
//SUBROUTINE RCARRY(RVEC,LENV)
void rcarry_(float *, int *);


double energy(int L,int S[L][L],int pbc[L][2]){
	double ene;
	int i,j;

	ene=0.;
	for(i=0;i<L;i++)
	{
		for(j=0;j<L;j++)
		{
			ene -= S[i][j]*S[pbc[i][1]][j]+S[i][pbc[j][1]];
		}
	}
	return ene;
}

int main(void){
	// ----- VARIABLE DEFINITION ---------------
	int L,i,j,k,irand;
	printf("L=");
	scanf("%d",&L);
	int nrand=L*L*3+24,seed;
	printf("Seed=");
	scanf("%d",&seed);
	float rrand[nrand];
	int pbc[L][2],S[L][L];
	double ene;
	// ------------------------------------------

	// ----- SPIN MATRIX DEFINITION -------------

	rcarin_(&seed,rrand,&nrand);
	rcarry_(rrand,&nrand);

	irand=0;
	for(i=0;i<L;i++)
	{
		for(j=0;j<L;j++)
		{
			if(rrand[nrand]<0.5)S[i][j]=1;
			else S[i][j]=-1;
			irand+=1;
		}
	}
	//------------------------------------------

	// ----- PBC DEFINITION --------------------
	pbc[0][0]=L-1;
	pbc[0][1]=1;
	pbc[L-1][0]=L-2;
	pbc[L-1][1]=0;
	for(i=1;i<L-1;i++)
	{
		pbc[i][0]=i-1;
		pbc[i][1]=i+1;
	}

	ene=energy(L,S,pbc);
	printf("%lf\n",ene);


	// ----- MAIN LOOP BEGINS ------------------
	rcarry_(rrand,&nrand);

	





	return 0;

}