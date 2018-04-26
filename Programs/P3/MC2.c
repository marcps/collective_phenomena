/*Marc Pascual - Universitat de Barcelona 2018 - Fenòmens Col·lectius i transicions de fase*/
/*###################################################################
  # name: MC2.c                                                     #
  # date: 20/4/18                                                   #
  #                                                                 #
  ###################################################################
  # programa que calcula promitjos de l'energia (E/N)(E²/N),i de la #
  # magnetització i les escriu en un arxiu de resultats.            #
  # La idea és que es corri el programa varies vegades per poder fer#
  # un gràfic amb gnuplot per diferents temperatures.               #
  #                                                                 #
  #       A l'hora de calcular els promitjos fem varis trucs:       #
  #		(1) No comptem varies dades inicials per donar temps#
  #                 a que la matriu d'spins s'estabilitzi           #
  #		(2) comptem cada 10 passes per evitar Correlacions  #
  #                 També tenim en compte el temps de CPU.          #
  ###################################################################*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
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
	//calculates the total energy of a spin system in 2-D
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
	//Correlacions
	int mctot,mcini,mcd;
	printf("[*]Nombre total de passes de MC: mctot=");
	scanf("%d",&mctot);
	printf("[*]Nombre de passes inicials que ens saltem: mcini=");
	scanf("%d",&mcini);
	printf("[*](10)mcd=");
	scanf("%d",&mcd);

	FILE *fp;
	clock_t start,end;
	double cpu_time_used;
	int L,i,j,k,p,irand,count;
	printf("[*]Longitud de la xarxa: L=");
	scanf("%d",&L);
	int nrand=mctot*3+24;
	float rrand[nrand];
	int pbc[L][2],S[L][L],itemp;
	double ene,magn,suma,de,vexp[4],temp,ftemp,epse,epsmcapv,capv_n,suscept,suscept_n;
	printf("[*]Temperatura inicial: temp=");
	scanf("%lf",&ftemp);

	//VARIABLES PROMITJOS:
	int sum=0;
	double sume=0.,sume2=0.,vare; //Promig de l'energia i l'energia al quadrat
	double summ=0.,summ2=0.,sumam=0.,varm; //Promig de m i m² i abs(m)

	//Promitjos sobre NLLAV llavors
	int nllav,nllav0,illav;
	printf("[*]Llavor inicial: nllav0=");
	scanf("%d",&nllav0);
	printf("[*]Nombre de llavors: nllav=");
	scanf("%d",&nllav);


	// ------------------------------------------
	fp=fopen("mc2-resultats.res","a");
	//CLOCK
	start = clock();
	for(itemp=0;itemp<=20;itemp++)
	{
		temp=ftemp+(double)itemp*ftemp/4.;
//############################# MAIN LOOP NLLAV BEGINS #########################################
		for(illav=nllav0;illav<nllav0+nllav;illav++)
		{
			//========= SPIN MATRIX DEFINITION ==========
			rcarin_(&illav,rrand,&nrand);//every time a new estocastic vector with new seed is ddefined
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
			//==========================================

			// ----- PBC DEFINITION --------------------
			//The first column refers to the actual position
			/*The second column refers to the previous position (0)
			Or to the next (1)*/
			pbc[0][0]=L-1;
			pbc[0][1]=1;
			pbc[L-1][0]=L-2;
			pbc[L-1][1]=0;
			for(i=1;i<L-1;i++)
			{
				//Trivial definition (but necessary)
				pbc[i][0]=i-1; //previous position
				pbc[i][1]=i+1; //Next position
			}

			ene=energy(L,S,pbc);
			printf("Energia inicial E=%lf\n",ene);
			//------------------------------------------
			rcarry_(rrand,&nrand);
			//Resetting counters
			irand=0; //counter for the estocastic vector "rrand"
			count=0; //counter for the average of physical prop (E,M)
		//=========--- MAIN LOOP BEGINS ----===================
			for(i=0;i<mctot;i++)
			{
				//This loop will test the change in energy of a change in spin
				//And reject it or not according to "de":

					//k,p are 2 randomly selected indexes of the spin matrix
					k=round((L-1)*rrand[irand]);
					irand++;
					p=round((L-1)*rrand[irand]);
					irand++;

					suma=S[k][pbc[p][0]]+S[k][pbc[p][1]]+
						S[pbc[k][0]][p]+S[pbc[k][1]][p];

					de=2.*S[k][p]*suma;

					if(de<0)S[k][p]=-S[k][p];
					else
					{
						if(rrand[irand]<exp(-de/temp))
						{
							S[k][p]=-S[k][p];
						}
						//If not, the value is not accepted
						irand++;
					}

					ene=energy(L,S,pbc);

					//-=========== PROMITJOS ===============-:
					if((count>mcini)&&(mcd*(count/mcd)==count)){
						magn=magne(L,S);

						sum++; //variable que compta quantes vegdes suma per fer el promig després
						sume += ene;
						sume2 += ene*ene;
						summ += magn;
						summ2 += magn*magn;
						sumam += fabs(magn);
					}
					//=======================================
					count++;
			}
		}
//######################### MAIN LOOP ENDS (NLLAV) ###################################################

		// Completing the average calculation:
		sume=sume/((double)sum);
		sume2=sume2/((double)sum);
		summ=summ/((double)sum);
		summ2=summ2/((double)sum);
		sumam=sumam/((double)sum);

		epse=(1./(L*L))*sqrt((sume2-sume*sume)/(double)sum);
		epsm=(1./(L*L))*sqrt((summ2-summ*summ)/sum);
		capv=((sume2-sume*sume)/temp)*((sume2-sume*sume)/temp);
		capv_n=capv/(L*L);
		suscept=((summ2-summ*summ)/temp)*((summ2-summ*summ)/temp);
		suscept_n=suscept/(L*L);

		vare=sume2-sume*sume;
		varm=summ2-summ*summ;

		end = clock();
		cpu_time_used = ((double)(end-start))/CLOCKS_PER_SEC;
		printf("\n----------------\n"
						"CPU TIME: %.1lf s\n"
						"-----------------\n",cpu_time_used);

		printf("\n\n\nNumber of operations=%d \nE/N=%.2lf\nE^2/N=%.2lf\n"
					"M/N= %.2lf\nM^2=%.2lf\nabs(M)/N=%.2lf\n\nVAR(E)=%.2lf\n"
					"VAR(M)=%.2lf\n\n",
					count,sume,sume2,summ,summ2,sumam,vare,varm);
		//Escriurem els resultats a l'arxiu per després fer una gràfica
		fprintf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf\n"
						,temp,sume,sume2,summ,summ2,sumam,vare,varm);
	}

	fclose(fp);
	return 0;
}
