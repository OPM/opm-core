/*===========================================================================
//
// File: call_petsc.h
//
// Created: 2014-05-07 10:21:21 CST
//
// Authors: Ming Liu  <miliu@statoil.com>
//==========================================================================*/
/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
  Copyright 2014 STATOIL ASA.

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
#ifndef OPM_CALL_PETSC_H_HEADER
#define OPM_CALL_PETSC_H_HEADER
#ifdef __cplusplus
extern "C" {
#endif
int 
call_Petsc(const int size, const int nonzeros, const int* ia, const int* ja, const double* sa, const double* b, double* x, int argc, char** argv, const int ksp_type, const int pc_type, const double rtol, const double atol, const double dtol, const int maxits, const int view_ksp);

#ifdef __cplusplus
}
#endif
#endif  /* OPM_CALL_PETSC_H_HEADER */
