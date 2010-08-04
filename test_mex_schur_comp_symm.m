run ../../startup
G = computeGeometry(cartGrid([2,1]));
rock.perm = ones(G.cells.num, 1);
BI = mex_ip_simple(G, rock);

nconn = diff(G.cells.facePos)
[S, r, F, L] = mex_schur_comp_symm(BI, nconn);

[i, j] = blockDiagIndex(nconn, nconn);

SS = sparse(double(G.cells.faces(i, 1)), ...
            double(G.cells.faces(j, 1)),  S );
R  = accumarray(G.cells.faces(:, 1), r);

x = SS \ R
