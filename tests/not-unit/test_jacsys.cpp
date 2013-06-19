#include "config.h"
#include <cstddef>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>

template <class Ostream, class Coll>
Ostream&
operator<<(Ostream& os, const Coll& c) {
	typedef typename Coll::value_type VT;

	os << "[ ";
	std::copy(c.begin(), c.end(), std::ostream_iterator<VT>(os, " "));
	os << "]";

	return os;
}

#if 0
int
main() {
	std::vector<int> a;

	a.push_back(0);
	a.push_back(1);

	std::cerr << "a[0] = " << a << '\n';

	for (::std::size_t i = 1; i < 7; ++i) {
		a.reserve(2 * a.size());
		a.insert(a.end(), a.begin(), a.end());

		std::cerr << "a[" << i << "] = " << a << '\n';
	}
}
#endif


#if 0
#include "VectorBlockAssembler.hpp"

using namespace Opm::LinAlgSupport;

int
main() {
	typedef std::vector<double>         stvd;

	typedef AccumulationNorm<SumAbs>    OneNorm;
	typedef AccumulationNorm<Euclid>    TwoNorm;
	typedef AccumulationNorm<MaxAbs>    MaxNorm;

	typedef NumericVector<stvd,OneNorm> VectorOne;
	typedef NumericVector<stvd,TwoNorm> VectorTwo;
	typedef NumericVector<stvd,MaxNorm> VectorInf;

	VectorOne v1;

	v1.resize(10, 1.0);

	std::cerr << v1 << '\n';
	std::cerr << "Norm = " << v1.norm() << '\n';

	VectorTwo v2;

	v2.resize(10, 1.0);

	std::cerr << v2 << '\n';
	std::cerr << "Norm = " << v2.norm() << '\n';

	VectorInf vi;

	vi.resize(10, 1.0);

	std::cerr << vi << '\n';
	std::cerr << "Norm = " << vi.norm() << '\n';
}
#endif

#if 0
#include <opm/core/transport/JacobianSystem.hpp>

//using namespace Opm::LinAlgSupport;
using namespace Opm::ImplicitTransportDefault;

int
main() {
	typedef std::vector<double>  STVD;

	NewtonVectorCollection<STVD> nvcoll;

	STVD blk(2, 1.0);

	nvcoll.setSize(blk.size(), 1);
	nvcoll.assembleBlock(blk.size(), 0, blk);

	std::cerr << "Residual  = " << nvcoll.residual () << '\n';
	std::cerr << "Increment = " << nvcoll.increment() << '\n';
	std::cerr << "Solution  = " << nvcoll.solution () << '\n';
}
#endif

#if 1
#include <opm/core/linalg/sparse_sys.h>

#include <opm/core/transport/JacobianSystem.hpp>
#include <opm/core/transport/CSRMatrixBlockAssembler.hpp>
#include <opm/core/transport/CSRMatrixUmfpackSolver.hpp>

using namespace Opm::ImplicitTransportDefault;

int
main() {
	typedef std::vector<double> STVD;
	typedef NewtonVectorCollection <STVD> NVecColl;

	JacobianSystem< struct CSRMatrix, NVecColl > JS;
	JS.setSize(1, 2, 2);
	std::cerr << "Residual = " << JS.vector().residual() << '\n';

	std::vector <int> conn(1);
	conn[0] = 0;
	JS.matasm().createBlockRow(0, conn, 1);

	double blk[] = { 1.0 };
	JS.matasm().assembleBlock(1, 0, 0, blk);
	JS.vector().assembleBlock(1, 0,    blk);

	conn[0] = 1;
	JS.matasm().createBlockRow(1, conn, 1);
	JS.matasm().assembleBlock(1, 1, 1, blk);
	JS.vector().assembleBlock(1, 1,    blk);

	const CSRMatrix& A = JS.matrix();
	csrmatrix_write_stream(&A, stderr);

	std::cerr << "Residual = " << JS.vector().residual() << '\n';

	Opm::ImplicitTransportLinAlgSupport::CSRMatrixUmfpackSolver linsolve;

	std::cerr << "Increment = " << JS.vector().increment() << '\n';
	linsolve.solve(A, JS.vector().residual(),
			JS.vector().writableIncrement());

	std::cerr << "Increment = " << JS.vector().increment() << '\n';

	std::cerr << "Solution = " << JS.vector().solution() << '\n';
	JS.vector().addIncrement();
	std::cerr << "Solution = " << JS.vector().solution() << '\n';
}
#endif
