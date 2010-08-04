function G = mcomputeGeometry(G)
   [fa,fc,fn,cc,cv] = mex_compute_geometry(G);
   G.faces.areas     = fa;
   G.faces.centroids = fc';
   G.faces.normals   = fn';
   G.cells.centroids = cc';
   G.cells.volumes   = cv;
end
