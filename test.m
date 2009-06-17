nx = 5;
ny = 7;
nz = 11;
%g=simpleGrdecl([nx, ny, nz], @(x) -0.055+0.11*x+0.011 );
g = makeModel3([200,220,30]);
%G=processGRDECL(g);
%clf,plotGrid(G);view(3);

g.ACTNUM=int32(g.ACTNUM);

grdecl = readGRDECL(fullfile(ROOTDIR, 'examples','grids','GSmodel.grdecl'));
