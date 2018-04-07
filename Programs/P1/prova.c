#include <stdio.h>
#include <stdlib.h>

//calls the Fortran SUBROUTINE RCARIN(IJKL,RVEC,LENV)
void rcarin_(int *,float *, int *);
//SUBROUTINE RCARRY(RVEC,LENV)
void rcarry_(float *, int *);

void writeconfig(int L, int **S){
  FILE *fp;
  int i,j;
  fp=fopen("P1-configuration.conf","w");
  for(i=0;i<L;i++)
  {
    for(j=0;j<L;j++)
    {
      if(S[i][j]==1)
      {
        fprintf(fp,"%d %d\n",i,j);
      }
    }
  }
  fclose(fp);
}

int main(int argc, char const *argv[]){
  int L=32,i,j,irand;
  int seed;
  printf("[*]Seed=");
  scanf("%d",&seed);
  //S as a pointer of pointers
  int **S;
  S=(int*)malloc(L*sizeof(int));
  for(i=0;i<L;i++)
  {
    S[i]=(int*)malloc(L*sizeof(int));
  }
  //-------------------------
  int nrand=L*L*3+24;
  float rrand[nrand];

  //Generation of the estocastic vector:
  rcarin_(&seed,rrand,&nrand);
  rcarry_(rrand,&nrand);

  //Generation of the spin matrix from the estocastic vector:
  irand=1; //COUNTER
  for(i=0;i<L;i++){
    for(j=0;j<L;j++){
      if(rrand[irand]<0.5){
        S[i][j]=1;
      }
      else S[i][j]=-1;
      irand += 1;
    }
  }
  writeconfig(L,S);
  return 0;
}
