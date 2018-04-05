#include <stdio.h>
#include <stdlib.h>

//calls the Fortran SUBROUTINE RCARIN(IJKL,RVEC,LENV)
void rcarin_(int *,float *, int *);
//SUBROUTINE RCARRY(RVEC,LENV)
void rcarry_(float *, int *);

int main(void){
  int L,i,j,seed,irand;
  printf("[*]Dimension of the spin MATRIX L=");
  scanf("%d",&L);
  printf("[*]Seed=");
  scanf("%d",&seed);
  //Number of estocastic numbers to be generated:
  int nrand=L*L*3+24;
  //Spin matrix
  int S[L][L];
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

  //Print out the matrix:
  for(i=0;i<L;i++){
    for(j=0;j<L;j++){
      if(S[i][j]==1)
      {
        printf(" %d ",S[i][j]);
      }
      else printf("%d ",S[i][j]);
    }
    printf("\n");
  }

}
