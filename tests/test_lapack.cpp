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
#include <vector>


int main()
{
    const int N = 5;
    const int nrhs = 1;
    double DU[N-1] = { 2.1,  -1.0,   1.9,   8.0 };
    double D[N]    = { 3.0,   2.3,  -5.0,  -0.9,   7.1 };
    double DL[N-1] = { 3.4,   3.6,   7.0,  -6.0 };
    double B[N] = { 2.7,  -0.5,   2.6,   0.6,   2.7 };
    // double B[N]    = { 2.7,  -0.5,   2.6,   0.6,   2.7 };
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

    //test of dgbsv_
    const int kl = 1;
    const int ku = 1;
    const int nrowAB = 2*kl + ku + 1;
    const int ldb = N;    
    const int ldab = nrowAB;
    int ipiv;
    std::vector<double> AB(nrowAB*N, 0.);
    std::vector<double> BB(N, 0.);
    double ABu[N] = { 0., 2.1,  -1.0,   1.9,   8.0};
    double ABd[N] = { 3.0,   2.3,  -5.0,  -0.9,   7.1 };
    double ABl[N] = { 3.4,   3.6,   7.0, - 6.0, 0. };
    for (int i = 0; i<N; ++i) {
        AB[kl + i*nrowAB] = ABu[i];
        AB[kl + 1 + i*nrowAB] = ABd[i];
        AB[kl + 2 + i*nrowAB] = ABl[i];
    }
    BB[0] = 2.7;
    BB[1] = -0.5;
    BB[2] = 2.6;
    BB[3] = 0.6;
    BB[4] = 2.7;
    std::vector<double>::iterator it;
    
    std::cout << "myvector contains:";
    for ( it=AB.begin() ; it < AB.end(); it++ ) {
        std::cout << " " << *it;
    }
    
    std::cout << std::endl;

    dgbsv_(&N, &kl, &ku, &nrhs, &AB[0], &ldab, &ipiv, &BB[0], &ldb, &info);
    if (info == 0) {
	for (int i = 0; i < N; ++i) {
	    std::cout << BB[i] << ' ';
	}
	std::cout << std::endl;
    } else {
	std::cerr << "Something went wrong in dgbsv_()\n";
    }
}

