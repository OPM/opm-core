run ../../startup
G = computeGeometry(cartGrid([200, 1]));
rock.perm = ones(G.cells.num, 1);
BI = mex_ip_simple(G, rock);

nconn = diff(G.cells.facePos);
conn  = G.cells.faces(:, 1);
[S, r, F, L] = mex_schur_comp_symm(BI, nconn, conn);

[i, j] = blockDiagIndex(nconn, nconn);

SS = sparse(double(conn(i)), double(conn(j)), S);
R  = accumarray(G.cells.faces(:, 1), r);

SS(1) = SS(1) * 2;
R([1, G.cells.num+1]) = [1, -1];

x = SS \ R;
[v, p] = mex_compute_press_flux(BI, x, nconn, conn, F, L);

plotCellData(G, p);
