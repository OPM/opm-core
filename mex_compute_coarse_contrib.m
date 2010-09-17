function varargout = mex_generate_coarse_contrib(varargin)
%Generate contributions to coarse-scale linsys using compiled C code.
%
% SYNOPSIS:
%   [cell_ip, Binv] = ...
%	mex_generate_coarse_contrib(BIf, Psi, p, pconn, dof_pos, blk_ncf)
%
% PARAMETERS:
%   BIf     - Fine-scale (inverse) inner product matrices.
%
%   Psi     - Basis function values.  Non-zeros only.  Ordered by blocks.
%
%   p       - Partition vector.
%
%   pconn   - Fine-scale indirection array into connection table.
%             Typically corresponds to G.cells.facePos.
%
%   dof_pos - Coarse-scale indirection array into (coarse) connection
%             table.  Roughly equivalent to CG.cells.facePos, but only
%             for active faces.
%
%   blk_ncf - Number of (fine-scale) half-contacts per block.
%
% RETURNS:
%   cell_ip - Cell contributions to coarse-scale (block) inner products.
%
%   Binv    - Coarse-scale inverse inner product.  Semantically
%             equivalent to 'BI' output of 'mex_ip_simple', but for the
%             coarse grid represented by 'p'.
%
% NOTE:
%   This function is only intended for testing/developing a compiled
%   language implementation of the MsMFE method.
%
% SEE ALSO:
%   mex_ip_simple.

%{
#COPYRIGHT#
%}

% $Date$
% $Revision$

   buildmex CFLAGS="\$CFLAGS -Wall -Wextra -ansi -pedantic           ...
        -Wformat-nonliteral -Wcast-align -Wpointer-arith             ...
        -Wbad-function-cast -Wmissing-prototypes -Wstrict-prototypes ...
        -Wmissing-declarations -Winline -Wundef -Wnested-externs     ...
        -Wcast-qual -Wshadow -Wconversion -Wwrite-strings            ...
        -Wno-conversion -Wchar-subscripts -Wredundant-decls"         ...
   ...
        -O -largeArrayDims -DTIME_LOCAL=1    ...
   ...
	mex_compute_coarse_contrib.c coarse_sys.c partition.c dfs.c ...
   ...
	-lmwlapack -lmwblas

   % Call MEX'ed edition.
   [varargout{1:nargout}] = mex_generate_coarse_contrib(varargin{:});
end
