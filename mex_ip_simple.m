function varargout = mex_ip_simple(varargin)

   buildmex -O -v -largeArrayDims ...
            mex_ip_simple.c mimetic.c                  ...
            -lmwlapack -lmwblas

   [varargout{1:nargout}] = mex_ip_simple(varargin{:});
end
