
  g=simpleGrdecl([2, 1, 4], @(x) -0.05+0.1*x+0.01 );
G=processGRDECL(g);
%clf,plotGrid(G);view(3);

g.ACTNUM=int32(g.ACTNUM);
