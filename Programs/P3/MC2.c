/*Marc Pascual - Universitat de Barcelona 2018 - Fenòmens Col·lectius i transicions de fase*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>


//calls the Fortran SUBROUTINE RCARIN(IJKL,RVEC,LENV)
//rcarin_(&seed,rrand,&nrand);
//rcarry_(rrand,&nrand);
void rcarin_(int *,float *, int *);
//SUBROUTINE RCARRY(RVEC,LENV)
void rcarry_(float *, int *);

double magne(int L, int S[L][L]){
	int i,j;
	double mag=0.0;
	for(i=0;i<L;i++)
	{
		for(j=0;j<L;j++)
		{
			mag+=S[i][j];
		}
	}
	return mag;
}

double energy(int L,int S[L][L],int pbc[L][2]){
	double ene;
	int i,j;

	ene=0.;
	for(i=0;i<L;i++)
	{
		for(j=0;j<L;j++)
		{
			ene += -S[i][j]*S[pbc[i][1]][j]
					-S[i][j]*S[i][pbc[j][1]];
		}
	}
	return ene;
}

int main(int argc, char const *argv[]){
	// ----- VARIABLE DEFINITION ---------------
	FILE *fp;
	int L,i,j,k,p,irand,count;
	printf("[*]L=");
	scanf("%d",&L);
	int nrand=L*L*3+24,seed;
	printf("[*]Seed=");
	scanf("%d",&seed);
	float rrand[nrand];
	int pbc[L][2],S[L][L];
	double ene,magn,suma,de,vexp[4],temp;
	printf("[*]Temperatura: temp=");
	scanf("%lf",&temp);

	//VARIABLES PROMITJOS:
	double sume=0.,sume2=0.,vare; //Promig de l'energia i l'energia al quadrat
	double summ=0.,summ2=0.,sumam=0.,varm; //Promig de m i m² i abs(m)


	// ------------------------------------------
	fp=fopen("Simul-temp2000-MCtot1000","w");

	// ----- SPIN MATRIX DEFINITION -------------
	rcarin_(&seed,rrand,&nrand);
	rcarry_(rrand,&nrand);

	irand=0;
	for(i=0;i<L;i++)
	{
		for(j=0;j<L;j++)
		{
			if(rrand[irand]<0.5)S[i][j]=1;
			else S[i][j]= -1;
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
	printf("Energia inicial E=%lf\n",ene);
	//------------------------------------------


	// ----- MAIN LOOP BEGINS ------------------
	rcarry_(rrand,&nrand);
	irand =0;
	count=0;
	for(i=0;i<L;i++)
	{
		for(j=0;j<L;j++)
		{
//This loop will test the change in energy of a change in spin
//And reject it or not according to de:

			//k i p són els índex triats aleaoriament de la matriu d'spins
			k= round(L*rrand[irand]);
			p=round(L*rrand[irand+1]);
			irand+=2;
			
			suma=S[k][pbc[p][0]]+S[k][pbc[p][1]]+
				S[pbc[k][0]][p]+S[pbc[k][1]][p];
			
			de=2.*S[k][p]*suma;

			if(de<0) S[k][p]=-S[k][p];
			else
			{
				if(rrand[irand]<exp(-de/temp))
				{
					S[k][p]=-S[k][p];
				}
				//If not, the value is not accepted
				irand+=1;
			}

			ene=energy(L,S,pbc);
			magn=magne(L,S);

			//######### PROMITJOS ##############:
			sume += ene;
			sume2 += ene*ene;
			summ += magn;
			summ2 += magn*magn;
			sumam += fabs(magn);
			//##################################
			count+=1;
		}
		
	}
	// Normalització dels promitjos:
	sume=sume/((double)count);
	sume2=sume2/((double)count);
	summ=summ/((double)count);
	summ2=summ2/((double)count);
	sumam=sumam/((double)count);

	vare=sume2-sume*sume;
	varm=summ2-summ*summ;
	printf("%d %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf\n",count,sume,sume2,summ,summ2,sumam,vare,varm );

	fclose(fp);
	return 0;

}