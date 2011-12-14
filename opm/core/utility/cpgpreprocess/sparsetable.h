/*===========================================================================
//
// File: sparsetable.h
//
// Created: Fri Jun 19 08:47:45 2009
//
// Author: Jostein R. Natvig <Jostein.R.Natvig@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================*/

/*
Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
Copyright 2009, 2010 Statoil ASA.

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

#ifndef OPENRS_SPARSETABLE_HEADER
#define OPENRS_SPARSETABLE_HEADER

typedef struct{
  int m;        /* number of rows          */
  int *ptr;     /* row pointer of size m+1 */
  int position; /* first position in ptr that is not filled. */

  int n;        /* size of data            */
  void *data;   /* sparse table data       */
  

} sparse_table_t;

void             free_sparse_table     (sparse_table_t *tab);
sparse_table_t   *malloc_sparse_table  (int m, int n, int datasz);
sparse_table_t   *realloc_sparse_table (sparse_table_t *tab, int m, int n, int datasz);

#endif /* OPENRS_SPARSETABLE_HEADER */
