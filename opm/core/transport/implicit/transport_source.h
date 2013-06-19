/*===========================================================================
//
// File: transport_source.h
//
// Created: 2011-10-05 19:58:53+0200
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

#ifndef OPM_TRANSPORT_SOURCE_H_HEADER
#define OPM_TRANSPORT_SOURCE_H_HEADER

#ifdef __cplusplus
extern "C" {
#endif

struct TransportSource {
    int nsrc;
    int cpty;

    int nphase;

    int    *cell;
    double *pressure;
    double *flux;
    double *saturation;
    double *surfvolume;
};


struct TransportSource *
create_transport_source(int nsrc, int nphase);

void
destroy_transport_source(struct TransportSource *src);

int
append_transport_source(int                     c,
                        int                     nphase,
                        double                  p,
                        double                  v,
                        const double           *sat,
                        const double           *z,
                        struct TransportSource *src);

void
clear_transport_source(struct TransportSource *src);

#ifdef __cplusplus
}
#endif

#endif  /* OPM_TRANSPORT_SOURCE_H_HEADER */
