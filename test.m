
g=simpleGrdecl([4, 2, 4], @(x) -0.055+0.11*x+0.011 );
G=processGRDECL(g);
%clf,plotGrid(G);view(3);

g.ACTNUM=int32(g.ACTNUM);
