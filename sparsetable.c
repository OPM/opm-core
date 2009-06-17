#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <mex.h>

#include "sparsetable.h"

void free_sparse_table (sparse_table_t *tab)
{
  if(tab->ptr)  free(tab->ptr);
  if(tab->data) free(tab->data);

  free(tab);
}

sparse_table_t *malloc_sparse_table(int m, int n, int datasz)
{
  sparse_table_t *tab = malloc(sizeof *tab);
  tab->m = m; 
  tab->n = n;
  tab->position = 0;
  if (!(tab->ptr  = malloc((m+1) * sizeof (*tab->ptr)))){
    fprintf(stderr, "Could not allocate space for sparse ptr\n");
    free_sparse_table(tab);
    return NULL;
  }
  
  
  if(!(tab->data = malloc(n * datasz))){
    fprintf(stderr, "Could not allocate space for sparse data(%d)\n", n);
    free_sparse_table(tab);
    return NULL;
  }    
  
  return tab;
}

sparse_table_t *realloc_sparse_table(sparse_table_t *tab, int m, int n, int datasz)
{  
  void *p = realloc(tab->ptr, (m+1) * sizeof (*tab->ptr));
  if (p){
    tab->ptr = p;
  }else{
    fprintf(stderr, "Could not reallocate space for sparse ptr\n");
    free_sparse_table(tab);
    return NULL;
  }
  
  p = realloc(tab->data, n * datasz);
  if(p){
    tab->data = p;
  }else{
    fprintf(stderr, "Could not reallocate space for sparse data(%d)\n", n);
    free_sparse_table(tab);
    return NULL;
  }
  
  tab->m = m; 
  tab->n = n;
  
  return tab;
}
