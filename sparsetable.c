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
  /* fprintf(stderr, "sizeof tab = %ld\n", sizeof (*tab->ptr)); */
  if (!(tab->ptr  = malloc((m+1) * sizeof (*tab->ptr)))){
    fprintf(stderr, "Could not allocate space for sparse ptr\n");
    free_sparse_table(tab);
    return NULL;
  }
  
  
  if(!(tab->data = malloc(n * datasz))){
    fprintf(stderr, "Could not allocate space for sparse data\n");
    free_sparse_table(tab);
    return NULL;
  }    
  
  return tab;
}

sparse_table_t *realloc_sparse_table(sparse_table_t *tab, int m, int n, int datasz)
{
  tab->m = m; 
  tab->n = n;

  if (!(tab->ptr  = realloc(tab->ptr, (m+1) * sizeof (*tab->ptr)))){
    fprintf(stderr, "Could not reallocate space for sparse ptr\n");
    free_sparse_table(tab);
    return NULL;
  }
  
  if(!(tab->data = realloc(tab->data, n * datasz))){
    fprintf(stderr, "Could not reallocate space for sparse data\n");
    free_sparse_table(tab);
    return NULL;
  }

  return tab;
}
