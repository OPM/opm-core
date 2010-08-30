function varargout = mex_generate_coarsegrid(varargin)
%Build coarse grid data structure from partition of existing grid (MEX)
%
% SYNOPSIS:
%   CG = mex_generate_coarsegrid(G, p)
%
% PARAMETERS:
%   G - Grid data structure as described in 'grid_structure'.
%
%   p - Partition vector as defined by, e.g., functions 'partitionUI' or
%       'partitionNonUniform'.
%
% RETURNS:
%   CG - Coarse grid data structure as described in 'generateCoarseGrid'.
%        There is, however, a number of subtle differences in the details
%        of this structure as compared to the pure M implementation.
%        Specifically, the coarse faces are numbered differently, and the
%        MEX implementation does not distinguish cardinal directions whence
%        the 'subFaces' function does not produce meaningful results on
%        outer faces.
%
% SEE ALSO:
%   generateCoarseGrid.

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
        -O -largeArrayDims ...
   ...
        mex_generate_coarsegrid.c coarse_conn.c mrst_abi.c

   % Call MEX'ed edition.
   [varargout{1:nargout}] = mex_generate_coarsegrid(varargin{:});
end
