#include <stdio.h>
#include <stdlib.h>

//calls the Fortran SUBROUTINE RCARIN(IJKL,RVEC,LENV)
void rcarin_(int *,float *, int *);
//SUBROUTINE RCARRY(RVEC,LENV)
void rcarry_(float *, int *);
