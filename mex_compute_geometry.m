function varargout = mex_compute_geometry(varargin)
   buildmex -O -largeArrayDims -lm ...
            mex_compute_geometry.c mimetic_geometry.c mrst_api.c 
            

   % Call MEX'ed edition.
   [varargout{1:nargout}] = mex_compute_geometry(varargin{:});
end
