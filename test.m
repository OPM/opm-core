nx = 6;
ny = 6;
nz = 11;
%g=simpleGrdecl([nx, ny, nz], @(x) 0.05+0.11*x+0.011 );
g = makeModel3([5,5,300]);
%G=processGRDECL(g);
%clf,plotGrid(G);view(3);

g.ACTNUM=int32(g.ACTNUM);

grdecl = readGRDECL(fullfile(ROOTDIR, 'examples','grids','GSmodel.grdecl'));
