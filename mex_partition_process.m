function varargout = mex_partition_process(varargin)
%Split disconnected blocks into new blocks using compiled C code.
%
% SYNOPSIS:
%   p = mex_partition_process(p, Neighbours)
%
% PARAMETERS:
%   p - Partition vector.  Should not contain any empty blocks.  Use
%       function 'mex_partition_compress' to remove empty blocks/bins.
%
%   Neighbours -
%        Neighbour definition.  An m-by-two numeric array of cell-to-cell
%        connections.  Often equal to the 'faces.neighbors' array of an MRST
%        'grid_structure', but general definitions are supported.
%
% RETURNS:
%   p - Updated partition vector.  Disconnected blocks of original partition
%       vector are assigned new block numbers.  Specifically, if block 'b'
%       contains 'ncomp' separate connected components, then ncomp-1 new
%       blocks are created from block 'b'.
%
% NOTE:
%   This implementation corresponds to setting Reconnect=FALSE in the
%   'processPartition' function when 'Neighbour' is G.faces.neighbor of a
%   grid_structure, 'G'.
%
% SEE ALSO:
%   processPartition.

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
        mex_partition_process.c partition.c dfs.c

   % Call MEX'ed edition.
   [varargout{1:nargout}] = mex_partition_process(varargin{:});
end
