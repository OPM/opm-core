function varargout = mex_partition_ui(varargin)
%Partition grid uniformly in logical space using compiled C code.
%
% SYNOPSIS:
%   p = mex_partition_ui(ix, fineDim, coarseDim)
%
% PARAMETERS:
%   ix        - Indices, e.g., G.cells.indexMap.
%
%   fineDim   - Cartesian dimensions of underlying index space (e.g.
%               G.cartDims)
%
%   coarseDim - Cartesian dimensions of requested index space.  Typically,
%               all(coarseDim <= fineDim).
%
% RETURNS:
%   p - Partition vector of SIZE(ix).
%
% SEE ALSO:
%   partitionUI.

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
        -O -largeArrayDims -DCOMPILING_FOR_MATLAB=1    ...
   ...
        mex_partition_ui.c partition.c

   % Call MEX'ed edition.
   [varargout{1:nargout}] = mex_partition_ui(varargin{:});
end
