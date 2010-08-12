run ../../startup
G = computeGeometry(cartGrid([200, 1], [1, 1]));
rock.perm = ones(G.cells.num, 1);

[BI, connPos, conns] = mex_ip_simple(G, rock);

[S, r, F, L] = mex_schur_comp_symm(BI, connPos, conns);

nconn  = diff(connPos);
[i, j] = blockDiagIndex(nconn, nconn);

SS = sparse(double(conns(i)), double(conns(j)), S);
R  = accumarray(conns, r);

SS(1) = SS(1) * 2;
R([1, G.cells.num+1]) = [1, -1];

x = SS \ R;
[v, p] = mex_compute_press_flux(BI, x, connPos, conns, F, L);

plotCellData(G, p);
