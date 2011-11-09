function varargout = processgrid(varargin)
%Compute grid topology and geometry from pillar grid description.
%
% SYNOPSIS:
%   G = processgrid(grdecl)
%   G = processgrid(grdecl,ztol)
%
% PARAMETERS:
%   grdecl - Raw pillar grid structure, as defined by function
%            'readGRDECL', with fields COORDS, ZCORN and, possibly, ACTNUM.
%   ztol   - tolerance for unique points 
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

[varargout{1:nargout}] = processgrid_mex(varargin{:});
varargout{1}.griddim = 3;
