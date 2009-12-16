function varargout = processgrid(varargin)
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
buildmex CFLAGS='$CFLAGS -Wall -fPIC' processgrid.c preprocess.c ...
         uniquepoints.c facetopology.c sparsetable.c mxgrdecl.c

% Call MEX edition.
[varargout{1:nargout}] = processgrid(varargin{:});
