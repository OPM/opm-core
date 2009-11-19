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

% Copyright 2009 SINTEF ICT, Applied Mathematics.
% Mex gateway by Jostein R. Natvig, SINTEF ICT.

% Build
disp('building processgrid.mex*')
op = pwd;
p  = which('processgrid');
ix = strfind(p, filesep);
p  = p(1:ix(end));
cd(p);
mex processgrid.c preprocess.c uniquepoints.c facetopology.c ...
    sparsetable.c mxgrdecl.c CFLAGS='$CFLAGS -Wall -fPIC'
cd(op)



%run
[varargout{1:nargout}] = processgrid(varargin{:});
