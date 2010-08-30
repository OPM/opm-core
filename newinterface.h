/*===========================================================================
//
// File: preprocess.h
//
// Created: Fri Jun 19 08:43:04 2009
//
// Author: Jostein R. Natvig <Jostein.R.Natvig@sintef.no>
//
// $Date: 2010-08-27 19:12:16 +0200 (Fri, 27 Aug 2010) $
//
// $Revision: 930 $
//
//==========================================================================*/

/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.
*/

#ifndef NEWINTERFACE_H
#define NEWINTERFACE_H
#include "preprocess.h"
#include "../reorder-utils/grid.h"

#ifdef __cplusplus
extern "C" {
#endif

   typedef struct {
      /* 
       *   Common grid definitions 
       */
      GRID_TOPOLOGY
      GRID_GEOMETRY
      
      /* 
       *   Special cornerpoint definitions 
       */
      int            cartdims[3];
      enum face_tag *face_tag;
      int            number_of_nodes_on_pillars;
      int           *cartesian_cell_index;
      
   } cornerpoint_grid_t;
   
   
   void preprocess         (const struct grdecl   *in, 
                            double                tol, 
                            cornerpoint_grid_t    *out);
   
   void free_cornerpoint_grid(cornerpoint_grid_t    *g);
   
#ifdef __cplusplus
}
#endif

#endif
