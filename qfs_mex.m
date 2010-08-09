cartDims = [500, 500];
physDims = [1000, 1000];
g = computeGeometry(cartGrid(cartDims, physDims));

rock.perm = ones(g.cells.num, 1);

nconn = diff(g.cells.facePos);
nconn([1,end]) = nconn([1, end]) + 1;
conn = zeros(sum(nconn), 1);

pos = cumsum([1; double(nconn)]);
ii  = mcolon(pos(1 : end-1), ...
             pos(1 : end-1) - 1 + double(diff(g.cells.facePos)));

conn(ii) = g.cells.faces(:,1);
conn([nconn(1), end]) = g.faces.num + (1:2);

BI = mex_ip_simple(g, rock, nconn, conn);
BI([nconn(1)^2, end]) = 1;   % Fake production indices...

[S, r, F, L] = mex_schur_comp_symm(BI, nconn, conn);

[i, j] = blockDiagIndex(nconn, nconn);
SS = sparse(double(conn(i)), double(conn(j)), S);
R  = accumarray(conn, r);

R(end-1 : end) = [1, -1];

lam = SS \ R;

[flux, press] = mex_compute_press_flux(BI, lam, nconn, conn, F, L);

plotCellData(g, press);
