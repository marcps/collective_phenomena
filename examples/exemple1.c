#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*Aquest exercici s'ha dut a terme amb l'exemple següent:
                p0(X) = (0,0.5,0.5,0)
                W( X->Y ) =(0,   0.5 ,0.333,  0)
                           (0.5,  0,  0.333,  0)
                           (0.5, 0.5,  0,     1)
                           (0,     0,  0.333, 0)
Assegurant que la matriu és estacionaria*/
int main(int argc, char const *argv[]) {
  int d;
  printf("Nombre de dimensions. d=");
  scanf("%d",&d);
  double v[d],g[d],p[d], W[d][d],temp;
  int i,j,k,N; //Nombre de passes
  FILE *fp;

  N=1;
  //INPUT
  for (i=0;i<d;i++){
    printf("p[%d]=",i+1 );
    scanf("%lf",&p[i]);
  }
  for (i=0;i<d;i++){
    for(j=0;j<d;j++){
      printf("W[%d][%d]=",i+1,j+1 );
      scanf("%lf",&W[i][j]);
    }
  }

  fp=fopen("resultats-ex1.txt","w");
  while(N<100){

    for (i=0;i<d;i++){
      g[i]=p[i];
    }
    //------------------------------
    //Multiplicació de la matriu pel vector k vegades
    for (k=0;k<N;k++){
      //Ho posem tot a zero
      temp=0.0;
      //Multiplicació de matrius
      for(i=0;i<d;i++){
        temp=0.0;
        for(j=0;j<d;j++){
          temp += g[j]*W[i][j];
        }
        v[i]=temp;
      }

      for (i=0;i<d;i++) {
        g[i]=v[i];
      }
    }
    //------------------------------

    //Escriurem els resultats a l'arxiu
    fprintf(fp,"%d ",N);
    for(i=0;i<d;i++)
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
