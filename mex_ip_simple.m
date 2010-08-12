function varargout = mex_ip_simple(varargin)
%Compute 'ip_simple' inner product values using compiled C code.
%
% SYNOPSIS:
%   [BI, connPos, conns] = mex_ip_simple(G, rock)
%   [BI, connPos, conns] = mex_ip_simple(G, rock, W)
%
% PARAMETERS:
%   G    - Grid data structure.
%
%   rock - Rock data structure.  Must contain a valid field 'perm'.
%
%   W    - Well data structure as defined by function 'addWell'.
%          OPTIONAL.  Only processed if present and valid.
%
% RETURNS:
%   BI      - An array of inner product values, ordered by the cells of the
%             input grid.  The array contains SUM(DIFF(connPos) .^ 2)
%             elements.
%
%   connPos - Indirection map of size [G.cells.num,1] into 'conns' table
%             (i.e., the connections or DOFs).  Specifically, the DOFs
%             connected to cell 'i' are found in the submatrix
%
%                  conns(connPos(i) : connPos(i + 1) - 1)
%
%             In the absence of wells, connPos == G.cells.facePos.
%
%   conns   - A (connPos(end)-1)-by-1 array of cell connections
%             (local-to-global DOF mapping in FEM parlance).  In the
%             absence of wells, conns == G.cells.faces(:,1).
%
% NOTE:
%   As the return value 'BI' is but a simple data array value, it must be
%   subsequently assembled into the 'S.BI' sparse matrix before being used
%   to solve a flow problem using, e.g., the 'solveIncompFlow' function.
%
%   Moreover, the 'solveIncompFlow' function expects its 'S' parameter to
%   specify a 'type' field which is consistent with the kind of matrix
%   stored within 'S'.  In the case of 'ip_simple', the 'type' must be the
%   string value 'hybrid'.
%
% EXAMPLE:
%   G = computeGeometry(processGRDECL(makeModel3([100, 60, 15])));
%   K = logNormLayers(G.cartDims, [10, 300, 40, 0.1, 100]);
%   rock.perm = bsxfun(@times, [1, 100, 0.1], K(:));
%   rock.perm = convertFrom(rock.perm(G.cells.indexMap, :), ...
%                           milli*darcy);
%
%   t0 = tic;
%   [BI, connPos, conns] = mex_ip_simple(G, rock);
%   toc(t0)
%
%   nconn = diff(connPos);
%   [i, j] = blockDiagIndex(nconn, nconn);
%
%   S = struct('BI', sparse(i, j, BI), 'type', 'hybrid', 'ip', 'ip_simple')
%
%   t0 = tic;
%   S2 = computeMimeticIP(G, rock)
%   toc(t0)
%
%   norm(S.BI - S2.BI, inf) / norm(S2.BI, inf)
%
% SEE ALSO:
%   computeMimeticIP, solveIncompFlow, blockDiagIndex.

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
        -O -largeArrayDims -DCOMPILING_FOR_MATLAB=1 ...
        mex_ip_simple.c mimetic.c mrst_api.c        ...
        -lmwlapack -lmwblas

   % Call MEX'ed edition.
   [varargout{1:nargout}] = mex_ip_simple(varargin{:});
end
