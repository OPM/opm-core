function varargout = mex_compute_press_flux(varargin)
%Derive pressure and flux from hybrid system using compiled C code.
%
% SYNOPSIS:
%   [v, p] = mex_compute_press_flux(BI, lam, nconn, conn, F, L)
%
% PARAMETERS:
%   BI    - Inner product values.
%
%   nconn - diff(G.cells.facePos)
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
%   BI = mex_ip_simple(G, rock);
%   toc(t0)
%
%   nconn = diff(G.cells.facePos);
%
%   [S, r, F, L] = mex_schur_comp(BI, nconn);
%
%   [i, j] = blockDiagIndex(diff(G.cells.facePos), ...
%                           diff(G.cells.facePos));
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

   buildmex -O -largeArrayDims -DCOMPILING_FOR_MATLAB=1 ...
            mex_compute_press_flux.c hybsys.c           ...
            -lmwlapack -lmwblas

   % Call MEX'ed edition.
   [varargout{1:nargout}] = mex_compute_press_flux(varargin{:});
end
