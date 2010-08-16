function varargout = mex_partition_compress(varargin)
%Partition grid uniformly in logical space using compiled C code.
%
% SYNOPSIS:
%   p = mex_partition_compress(p)
%
% PARAMETERS:
%   p - Original partition vector, may contain empty coarse blocks.
%
% RETURNS:
%   p - Updated partition vector.
%       Renumbered so as to remove empty coarse blocks.
%
% SEE ALSO:
%   compressPartition, partitionUI, mex_partition_ui.
%
% EXAMPLE:
%   p = partitionCartGrid([4, 4, 1], [2, 2, 2]);
%   [p, compressPartition(p), mex_partition_compress(p)]

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
        mex_partition_compress.c partition.c

   % Call MEX'ed edition.
   [varargout{1:nargout}] = mex_partition_compress(varargin{:});
end
