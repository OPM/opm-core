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

namespace {
    struct BandMatrixCoeff
    {
        BandMatrixCoeff(int N, int ku, int kl) : ku_(ku), kl_(kl), nrow_(2*kl + ku + 1), N_(N) {
        }


        // compute the position where to store the coefficient of a matrix A_{i,j} (i,j=0,...,N-1)
        // in a array which is sent to the band matrix solver of LAPACK.

        int operator ()(int i, int j) {
            return kl_ + ku_ + i - j + j*nrow_;
        }

        const int ku_;
        const int kl_;
        const int nrow_;
        const int N_;
    };
} //end anonymous namespace

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
    int kl = 1;
    int ku = 1;
    int nrowAB = 2*kl + ku + 1;
    int ldb = N;
    int ldab = nrowAB;
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


    dgbsv_(&N, &kl, &ku, &nrhs, &AB[0], &ldab, &ipiv, &BB[0], &ldb, &info);
    if (info == 0) {
	for (int i = 0; i < N; ++i) {
	    std::cout << BB[i] << ' ';
	}
	std::cout << std::endl;
    } else {
	std::cerr << "Something went wrong in dgbsv_()\n";
    }

    int NN = 3;
    kl = 1;
    ku = 1;
    nrowAB = 2*kl + ku + 1;
    ldb = NN;
    ldab = nrowAB;
    AB.assign(NN*nrowAB, 0.);
    BB.assign(NN, 0.);
    BandMatrixCoeff bmc(NN, ku, kl);
    AB[bmc(0, 0)] = 1.;
    AB[bmc(1, 0)] = 1.;
    AB[bmc(0, 1)] = 1.;
    AB[bmc(1, 1)] = 2.;
    AB[bmc(2, 1)] = 2.;
    AB[bmc(1, 2)] = 1.;
    AB[bmc(2, 2)] = 4.;
    BB[0] = 3.;
    BB[1] = 8.;
    BB[2] = 16.;

    // std::vector<double>::iterator it;
    // std::cout << "myvector contains:";
    // for ( it=AB.begin() ; it < AB.end(); it++ ) {
    //     std::cout << " " << *it;
    // }
    // std::cout << std::endl;

    dgbsv_(&NN ,&kl, &ku, &nrhs, &AB[0], &ldab, &ipiv, &BB[0], &ldb, &info);
    if (info == 0) {
        for (int i = 0; i < NN; ++i) {
            std::cout << BB[i] << ' ';
        }
        std::cout << std::endl;
    } else {
        std::cerr << "Something went wrong in dgbsv_()\n";
    }

}

