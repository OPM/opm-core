/*===========================================================================
//
// File: sparsetable.c
//
// Created: Fri Jun 19 08:48:04 2009
//
// Author: Jostein R. Natvig <Jostein.R.Natvig@sintef.no>
//
// $Date$
//
// $Revision$
//
//==========================================================================*/

/*
Copyright 2009 SINTEF ICT, Applied Mathematics.
Copyright 2009 Statoil ASA.

This file is part of The Open Reservoir Simulator Project (OpenRS).

OpenRS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or 
(at your option) any later version.

OpenRS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with OpenRS.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>


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
