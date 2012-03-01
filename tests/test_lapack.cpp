/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

  This file is part of the Open Porous Media project (OPM).

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

#include <opm/core/linalg/blas_lapack.h>
#include <iostream>

int main()
{
    const int N = 5;
    const int nrhs = 1;
    double DU[N-1] = { 2.1,  -1.0,   1.9,   8.0 };
    double D[N]    = { 3.0,   2.3,  -5.0,  -0.9,   7.1 };
    double DL[N-1] = { 3.4,   3.6,   7.0,  -6.0 };
    double B[N]    = { 2.7,  -0.5,   2.6,   0.6,   2.7 };
    int info = 0;
    dgtsv_(&N, &nrhs, DL, D, DU, B, &N, &info);
    if (info == 0) {
	for (int i = 0; i < N; ++i) {
	    std::cout << B[i] << ' ';
	}
	std::cout << std::endl;
    } else {
	std::cerr << "Something went wrong in dgtsv_()\n";
    }
}
