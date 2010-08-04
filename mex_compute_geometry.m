function varargout = mex_compute_geometry(varargin)
   buildmex -O -largeArrayDims                           ...
            mex_compute_geometry.c geometry.c mrst_api.c ...
            -lm

   % Call MEX'ed edition.
   [varargout{1:nargout}] = mex_compute_geometry(varargin{:});
end
