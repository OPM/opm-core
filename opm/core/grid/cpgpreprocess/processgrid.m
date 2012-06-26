function G = processgrid(varargin)
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

   G = processgrid_mex(varargin{:});
   G.griddim = 3;
   G = splitDisconnectedGrid(G, false);   
end

function G = splitDisconnectedGrid(G, verbose)
   % Check if grid is connected
   ix = all(G.faces.neighbors~=0, 2);
   I  = [G.faces.neighbors(ix,1);G.faces.neighbors(ix,2)];
   J  = [G.faces.neighbors(ix,2);G.faces.neighbors(ix,1)];
   N  = double(max(G.faces.neighbors(:)));
   A  = sparse(double(I),double(J),1,N,N)+speye(N);
   clear ix I J
   [a,b,c,d]=dmperm(A);                                %#ok
   clear A b d
   if numel(c) > 2,
      dispif(verbose, '\nGrid has %d disconnected components\n', ...
         numel(c)-  1);
      % Partition grid into connected subgrids
      for i = 1:numel(c) - 1,
         g(i)  = extractSubgrid(G, a(c(i):c(i+1)-1));  %#ok
         sz(i) = g(i).cells.num;                       %#ok
         g(i).cartDims = G.cartDims;       %#ok
      end
      
      % Return largest (in number of cells) grid first
      [i,i] = sort(-sz);                                                 %#ok
      G     = g(i);
   end
end
