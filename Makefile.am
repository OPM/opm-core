ACLOCAL_AMFLAGS = -I m4

# Recurse to build examples after libraries
SUBDIRS = . tests examples

# Declare libraries
lib_LTLIBRARIES = libopmcore.la

libopmcore_la_SOURCES = 				\
opmcore/utility/cart_grid.c				\
opmcore/utility/cpgpreprocess/mxgrdecl.c		\
opmcore/utility/cpgpreprocess/geometry.c		\
opmcore/utility/cpgpreprocess/preprocess.c		\
opmcore/utility/cpgpreprocess/processgrid.c		\
opmcore/utility/cpgpreprocess/readvector.cpp		\
opmcore/utility/cpgpreprocess/cgridinterface.c		\
opmcore/utility/cpgpreprocess/sparsetable.c		\
opmcore/utility/cpgpreprocess/facetopology.c		\
opmcore/utility/cpgpreprocess/scaneclipsedeck.c		\
opmcore/utility/cpgpreprocess/uniquepoints.c		\
opmcore/linalg/sparse_sys.c				\
opmcore/linalg/call_umfpack.c				\
opmcore/pressure/cfsh.c					\
opmcore/pressure/flow_bc.c				\
opmcore/pressure/well.c					\
opmcore/pressure/fsh_common_impl.c			\
opmcore/pressure/fsh.c					\
opmcore/pressure/tpfa/ifs_tpfa.c			\
opmcore/pressure/tpfa/compr_source.c			\
opmcore/pressure/tpfa/cfs_tpfa.c			\
opmcore/pressure/tpfa/compr_bc.c			\
opmcore/pressure/tpfa/compr_quant.c			\
opmcore/pressure/tpfa/compr_quant_general.c		\
opmcore/pressure/tpfa/cfs_tpfa_residual.c		\
opmcore/pressure/tpfa/trans_tpfa.c			\
opmcore/pressure/msmfem/coarse_conn.c			\
opmcore/pressure/msmfem/partition.c			\
opmcore/pressure/msmfem/hash_set.c			\
opmcore/pressure/msmfem/ifsh_ms.c			\
opmcore/pressure/msmfem/dfs.c				\
opmcore/pressure/msmfem/coarse_sys.c			\
opmcore/pressure/ifsh.c					\
opmcore/pressure/mimetic/mimetic.c			\
opmcore/pressure/mimetic/hybsys_global.c		\
opmcore/pressure/mimetic/hybsys.c			\
opmcore/transport/spu_explicit.c			\
opmcore/transport/spu_implicit.c			\
opmcore/transport/transport_source.c

nobase_include_HEADERS =				\
opmcore/utility/cpgpreprocess/readvector.hpp		\
opmcore/utility/cpgpreprocess/uniquepoints.h		\
opmcore/utility/cpgpreprocess/preprocess.h		\
opmcore/utility/cpgpreprocess/sparsetable.h		\
opmcore/utility/cpgpreprocess/cgridinterface.h		\
opmcore/utility/cpgpreprocess/geometry.h		\
opmcore/utility/cpgpreprocess/facetopology.h		\
opmcore/utility/cpgpreprocess/grid.h			\
opmcore/utility/cpgpreprocess/grdecl.h			\
opmcore/utility/cart_grid.h				\
opmcore/well.h						\
opmcore/linalg/sparse_sys.h				\
opmcore/linalg/call_umfpack.h				\
opmcore/linalg/blas_lapack.h				\
opmcore/fluid/SimpleFluid2p.hpp				\
opmcore/grid.h						\
opmcore/pressure/fsh.h					\
opmcore/pressure/HybridPressureSolver.hpp		\
opmcore/pressure/TPFAPressureSolver.hpp			\
opmcore/pressure/tpfa/compr_quant.h			\
opmcore/pressure/tpfa/compr_source.h			\
opmcore/pressure/tpfa/compr_bc.h			\
opmcore/pressure/tpfa/cfs_tpfa.h			\
opmcore/pressure/tpfa/compr_quant_general.h		\
opmcore/pressure/tpfa/trans_tpfa.h			\
opmcore/pressure/tpfa/cfs_tpfa_residual.h		\
opmcore/pressure/tpfa/ifs_tpfa.h			\
opmcore/pressure/TPFACompressiblePressureSolver.hpp	\
opmcore/pressure/GridAdapter.hpp			\
opmcore/pressure/msmfem/coarse_sys.h			\
opmcore/pressure/msmfem/hash_set.h			\
opmcore/pressure/msmfem/ifsh_ms.h			\
opmcore/pressure/msmfem/dfs.h				\
opmcore/pressure/msmfem/partition.h			\
opmcore/pressure/msmfem/coarse_conn.h			\
opmcore/pressure/flow_bc.h				\
opmcore/pressure/fsh_common_impl.h			\
opmcore/pressure/mimetic/hybsys_global.h		\
opmcore/pressure/mimetic/hybsys.h			\
opmcore/pressure/mimetic/mimetic.h			\
opmcore/transport/ImplicitTransport.hpp			\
opmcore/transport/CSRMatrixUmfpackSolver.hpp		\
opmcore/transport/CSRMatrixBlockAssembler.hpp		\
opmcore/transport/SimpleFluid2pWrapper.hpp		\
opmcore/transport/SinglePointUpwindTwoPhase.hpp		\
opmcore/transport/NormSupport.hpp			\
opmcore/transport/spu_explicit.h			\
opmcore/transport/JacobianSystem.hpp			\
opmcore/transport/spu_implicit.h			\
opmcore/transport/ImplicitAssembly.hpp			\
opmcore/transport/transport_source.h
