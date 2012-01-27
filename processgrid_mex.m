function varargout = processgrid_mex(varargin)
%Compute grid topology and geometry from pillar grid description.
%
% SYNOPSIS:
%   G = processGRDECL(grdecl)
%
% PARAMETERS:
%   grdecl - Raw pillar grid structure, as defined by function
%            'readGRDECL', with fields COORDS, ZCORN and, possibly, ACTNUM.
%
% RETURNS:
%   G      - Valid grid definition containing connectivity, cell
%            geometry, face geometry and unique nodes.
%
% EXAMPLE:
%   G = processgrid(readGRDECL('small.grdecl'));
%   plotGrid(G); view(10,45);
%
% SEE ALSO:
%   processGRDECL, readGRDECL, deactivateZeroPoro, removeCells.

%{
#COPYRIGHT#
%}

% $Date$
% $Revision$

% Copyright 2009 SINTEF ICT, Applied Mathematics.
% Mex gateway by Jostein R. Natvig, SINTEF ICT.

% $Date$
% $Revision$

   % Build MEX edition of same.
   %

   v = version;v = v([1,3]);

   CFLAGS = {'CFLAGS="\$CFLAGS', '-g', '-Wall', '-Wextra', '-ansi', ...
             '-pedantic', '-Wformat-nonliteral',  '-Wcast-align', ...
             '-Wpointer-arith', '-Wbad-function-cast', ...
             '-Wmissing-prototypes ', '-Wstrict-prototypes', ...
             '-Wmissing-declarations', '-Winline', '-Wundef', ...
             '-Wnested-externs', '-Wcast-qual', '-Wshadow', ...
             '-Wconversion', '-Wwrite-strings', '-Wno-conversion', ...
             '-Wchar-subscripts', '-Wredundant-decls"'};

   SRC = {'processgrid.c', 'preprocess.c', 'uniquepoints.c', ...
          'facetopology.c', 'sparsetable.c', 'mxgrdecl.c'};
       
   INCLUDE = {};
   
   OPTS = {'-output', ['processgrid_mex.', mexext], ...
      '-largeArrayDims', ['-DMATLABVERSION=', v], '-g'};
   
   buildmex(CFLAGS{:}, INCLUDE{:}, OPTS{:}, SRC{:})

   % Call MEX edition.
   [varargout{1:nargout}] = processgrid_mex(varargin{:});
end

