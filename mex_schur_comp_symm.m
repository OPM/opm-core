function varargout = mex_schur_comp_symm(varargin)
%Compute hybrid system component matrices using compiled C code.
%
% SYNOPSIS:
%   [S, r, F, L] = mex_schur_comp_symm(BI, nconn)
%
% PARAMETERS:
%   BI    - Inner product values.
%
%   nconn - Number of connections per cell.  Often coincides with
%           DIFF(G.cells.facePos), but may be larger if any cells are
%           perforated by one or more wells.
%
% RETURNS:
%   S - A SUM(nconn .^ 2)-by-1 array of unassembled system matrix values,
%       ordered by cells.
%
%   r - A SUM(nconn)-by-1 array of unassemble system rhs values, ordered by
%       cells.
%
%   F - A SUM(nconn)-by-1 array of C'*inv(B) values, ordered by cells.
%
%   L - A G.cells.num-by-1 array of C'*inv(B)*C values, ordered by cells.
%
% EXAMPLE:
%   G = computeGeometry(processGRDECL(makeModel3([100, 60, 15])));
%   K = logNormLayers(G.cartDims, [10, 300, 40, 0.1, 100]);
%   rock.perm = bsxfun(@times, [1, 100, 0.1], K(:));
%   rock.perm = convertFrom(rock.perm(G.cells.indexMap, :), ...
%                           milli*darcy);
%
%   nconn = diff(G.cells.facePos);
%   conn  = G.cells.faces(:,1);
%
%   BI = mex_ip_simple(G, rock, nconn, conn);
%
%   t0 = tic;
%   [S, r, F, L] = mex_schur_comp_symm(BI, nconn);
%   toc(t0)
%
%   [i, j] = blockDiagIndex(nconn, nconn);
%
%   SS = sparse(double(conn(i)), double(conn(j)), S);
%   R  = accumarray(conn, r);
%
%   lam = SS \ R;
%
% SEE ALSO:
%   mex_ip_simple, mex_compute_press_flux.

%{
#COPYRIGHT#
%}

% $Date$
% $Revision$

   buildmex -O CFLAGS="\$CFLAGS -W -Wall -pedantic -W -Wformat-nonliteral ...
        -Wcast-align -Wpointer-arith -Wbad-function-cast ...
        -Wmissing-prototypes -Wstrict-prototypes ...
        -Wmissing-declarations -Winline -Wundef -Wnested-externs...
        -Wcast-qual -Wshadow -Wconversion -Wwrite-strings...
        -Wno-conversion -Wchar-subscripts -Wredundant-decls" -largeArrayDims -DCOMPILING_FOR_MATLAB=1    ...
            mex_schur_comp_symm.c hybsys.c call_umfpack.c  ...
            -lmwlapack -lmwblas -lmwumfpack -lmwamd -I/work/include/SuiteSparse/

   % Call MEX'ed edition.
   [varargout{1:nargout}] = mex_schur_comp_symm(varargin{:});
end
