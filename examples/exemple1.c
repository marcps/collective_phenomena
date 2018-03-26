#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main(int argc, char const *argv[]) {
  double v[4],g[4],p[4], W[4][4],temp;
  int i,j,k,N; //Nombre de passes
  FILE *fp;

  N=1;
  //INPUT
  for (i=0;i<4;i++){
    printf("p[%d]=",i+1 );
    scanf("%lf",&p[i]);
  }
  for (i=0;i<4;i++){
    for(j=0;j<4;j++){
      printf("W[%d][%d]=",i+1,j+1 );
      scanf("%lf",&W[i][j]);
    }
  }

  fp=fopen("resultats-ex1.txt","w");
  while(N<100){

    for (i=0;i<4;i++){
      g[i]=p[i];
    }
    //------------------------------
    //Multiplicació de la matriu pel vector k vegades
    for (k=0;k<N;k++){
      //Ho posem tot a zero
      temp=0.0;
      //Multiplicació de matrius
      for(i=0;i<4;i++){
        temp=0.0;
        for(j=0;j<4;j++){
          temp += g[j]*W[i][j];
        }
        v[i]=temp;
      }

      for (i=0;i<4;i++) {
        g[i]=v[i];
      }
    }
    //------------------------------

    //Escriurem els resultats a l'arxiu
    fprintf(fp,"%d ",N);
    for(i=0;i<4;i++)
    {
      fprintf(fp,"%.10lf ",g[i]);
    }
    fprintf(fp,"\n");
    //end - Escriptura a l'arxiu
    N += 1;
  }
  //END WHILE
  fclose(fp);
  return 0;
}
