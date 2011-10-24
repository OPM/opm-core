/*===========================================================================
//
// File: compr_bc.h
//
// Created: 2011-10-24 15:58:16+0200
//
// Authors: Ingeborg S. Ligaarden <Ingeborg.Ligaarden@sintef.no>
//          Jostein R. Natvig     <Jostein.R.Natvig@sintef.no>
//          Halvor M. Nilsen      <HalvorMoll.Nilsen@sintef.no>
//          Atgeirr F. Rasmussen  <atgeirr@sintef.no>
//          Bård Skaflestad       <Bard.Skaflestad@sintef.no>
//
//==========================================================================*/


/*
  Copyright 2011 SINTEF ICT, Applied Mathematics.
  Copyright 2011 Statoil ASA.

  This file is part of the Open Porous Media Project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_COMPR_BC_H_HEADER
#define OPM_COMPR_BC_H_HEADER

#ifdef __cplusplus
extern "C" {
#endif

enum compr_bc_type { BC_PRESSURE, BC_RESV };

struct compr_bc {
    int                 nbc;
    int                 cpty;

    int                 nphases;

    enum compr_bc_type *type;
    int                *face;

    double             *press;
    double             *flux;
    double             *saturation;
};


struct compr_bc *
compr_bc_allocate(int np, int nbc);

void
compr_bc_deallocate(struct compr_bc *bc);

int
compr_bc_append(enum compr_bc_type  type,
                int                 f   ,
                int                 np  ,
                double              p   ,
                double              v   ,
                const double       *sat ,
                struct compr_bc    *bc  );

void
compr_bc_clear(struct compr_bc *bc);


#ifdef __cplusplus
}
#endif

#endif  /* OPM_COMPR_BC_H_HEADER */
