function varargout = mex_compute_press_flux(varargin)
%Derive pressure and flux from hybrid system using compiled C code.
%
% SYNOPSIS:
%   [v, p] = mex_compute_press_flux(BI, lam, nconn, conn, F, L)
%
% PARAMETERS:
%   BI    - Inner product values.  Typically computed using function
%           'mex_ip_simple'.
%
%   lam   - Interface pressure values.  One scalar value for each face in
%           the discretised reservoir model.
%
%   nconn - Number of connections per cell.  Typically equals the result of
%
%               diff(G.cells.facePos)
%
%   conn  - Actual connections per cell.  Typically equals the result of
%
%               G.cells.faces(:,1)
%
%   F     - Second-to-last return value from 'mex_schur_comp_symm'.
%
%   L     - Last return value from 'mex_schur_comp_symm'.
%
% RETURNS:
%   v - A SUM(nconn)-by-1 array of half-contact fluxes, ordered by cells.
%
%   p - A NUMEL(nconn)-by-1 array of cell pressure values.
%
% NOTE:
%   This function is the MEX'ed equivalent to the post-processing of
%   function 'schurComplementSymm'.  Note furthermore that this function
%   can only be used in conjunction with function 'mex_schur_comp_symm'.
%
% EXAMPLE:
%   G = computeGeometry(processGRDECL(makeModel3([100, 60, 15])));
%   K = logNormLayers(G.cartDims, [10, 300, 40, 0.1, 100]);
%   rock.perm = bsxfun(@times, [1, 100, 0.1], K(:));
%   rock.perm = convertFrom(rock.perm(G.cells.indexMap, :), ...
%                           milli*darcy);
%
%   nconn = diff(G.cells.facePos);
%   conn  = double(G.cells.faces(:,1));
%
%   BI = mex_ip_simple(G, rock);
%
%   [i, j] = blockDiagIndex(nconn, nconn);
%   [S, r, F, L] = mex_schur_comp_symm(BI, nconn);
%
%   SS = sparse(conn(i), conn(j), S);
%   R = accumarray(conn, r);
%
%   lam = SS \ R;
%
%   t0 = tic;
%   [v, p] = mex_compute_press_flux(BI, lam, nconn, conn, F, L);
%   toc(t0)
%
% SEE ALSO:
%   mex_ip_simple, mex_schur_comp_symm.

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
        mex_compute_press_flux.c hybsys.c           ...
        -lmwlapack -lmwblas

   % Call MEX'ed edition.
   [varargout{1:nargout}] = mex_compute_press_flux(varargin{:});
end
