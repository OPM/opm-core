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

#ifndef CGRIDINTERFACE_H
#define CGRIDINTERFACE_H

#include <opm/core/grid.h>
#include <opm/core/utility/cpgpreprocess/preprocess.h>

#ifdef __cplusplus
extern "C" {
#endif

    struct UnstructuredGrid *
    preprocess (const struct grdecl *in, double tol);

    void compute_geometry     (struct UnstructuredGrid *g);
    
    void free_grid(struct UnstructuredGrid *g);
#ifdef __cplusplus
}
#endif

#endif  /* CGRIDINTERFACE_H */
