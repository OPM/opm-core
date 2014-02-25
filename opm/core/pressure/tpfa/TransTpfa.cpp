#include "TransTpfa.hpp"
#include <opm/core/grid.h>

/* Ecplicitly initialize UnstructuredGrid versions */

template void tpfa_htrans_compute(const UnstructuredGrid*, const double*, double*);

template void tpfa_trans_compute(const UnstructuredGrid*, const double*, double*);

template void tpfa_eff_trans_compute(const UnstructuredGrid*, const double*, const double*, double*);
