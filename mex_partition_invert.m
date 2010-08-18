function varargout = mex_partition_invert(varargin)
%Invert cell-to-block map (creating block-to-cell) using compiled C code.
%
% SYNOPSIS:
%   [pb2c, b2c]      = mex_partition_invert(p)
%   [pb2c, b2c, loc] = mex_partition_invert(p)
%
% PARAMETERS:
%   p - Partition vector.  Should not contain any empty blocks.  Use
%       function 'mex_partition_compress' to remove empty blocks/bins.
%
% RETURNS:
%   pb2c - Indirection map of size [MAX(p) + 1, 1] into the 'b2c' map
%          array.  Specifically, the cells of block 'b' are stored in
%
%                b2c(pb2c(b) : pb2c(b + 1) - 1)
%
%   b2c  - Inverse cell map.  The entries in pb2c(b):pb2c(b+1)-1 correspond
%          to the result of FIND(p == b).
%
%   loc  - Local index within a block/bin.  Specifically,
%
%             loc(i) == FIND(b2c == i) - pb2c(p(i)) + 1
%
%          OPTIONAL.  Only returned (and computed) if specifically
%          requested.
%
% SEE ALSO:
%   mex_partition_ui, mex_partition_compress.

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
        mex_partition_invert.c partition.c

   % Call MEX'ed edition.
   [varargout{1:nargout}] = mex_partition_invert(varargin{:});
end
