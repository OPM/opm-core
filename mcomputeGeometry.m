function G = mcomputeGeometry(G)
%Compute geometric primitves using compiled C code.
%
% SYNOPSIS:
%   G = mcomputeGeometry(G)
%
% PARAMETERS:
%   G - A grid_structure.
%
% RETURNS:
%   G - An updated grid_structure containing areas, volumes, normals and
%       centroids.
%
% SEE ALSO:
%   computeGeometry, grid_structure.

%{
#COPYRIGHT#
%}

% $Date$
% $Revision$

   if numel(G) > 0,
      if ~isfield(G(1), 'type'),
         warning(msgid('GridType:Unknown'),                         ...
                ['Input grid has no known type. ',                  ...
                 'I''ll assume it arose from the primordial soup...']);

         [ G(:).type ] = deal( {'Primordial Soup'} );
      end

      for k = 1 : numel(G),
         [fa,fc,fn,cc,cv] = mex_compute_geometry(G(k));
         G(k).faces.areas     = fa;
         G(k).faces.centroids = fc';
         G(k).faces.normals   = fn';
         G(k).cells.centroids = cc';
         G(k).cells.volumes   = cv;
         
         G(k).type = [ G(k).type, { mfilename } ];
      end
   end
end
