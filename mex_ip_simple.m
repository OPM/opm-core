function varargout = mex_ip_simple(varargin)
%Compute 'ip_simple' inner product values using compiled C code.
%
% SYNOPSIS:
%   BI = mex_ip_simple(G, rock)
%
% PARAMETERS:
%   G    - Grid data structure.
%
%   rock - Rock data structure.  Must contain a valid field 'perm'.
%
% RETURNS:
%   BI   - A SUM(DIFF(G.cells.facePos) .^ 2)-by-1 array of inner product
%          values, ordered by the cells of the input grid.
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

   buildmex -O -v -largeArrayDims -DCOMPILING_FOR_MATLAB=1 ...
            mex_ip_simple.c mimetic.c                      ...
            -lmwlapack -lmwblas

   % Call MEX'ed edition.
   [varargout{1:nargout}] = mex_ip_simple(varargin{:});
end
