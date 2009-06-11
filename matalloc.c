#include <stdlib.h>
#include "matalloc.h"

double **dmatalloc(int m, int n) 
{
  double **ret = malloc(m*sizeof *ret);			
  ret[0] = malloc(m*n*sizeof *ret[0]);			
  int i;
  for (i=1; i<m; ++i) {
    ret[i] = ret[i-1]+n*sizeof *ret[0];
  }
  return ret;
}

int **imatalloc(int m, int n) 
{
  int **ret = malloc(m*sizeof *ret);			
  ret[0] = malloc(m*n*sizeof *ret[0]);			
  int i;
  for (i=1; i<m; ++i) {
    ret[i] = ret[i-1]+n*sizeof *ret[0];
  }
  return ret;
}
