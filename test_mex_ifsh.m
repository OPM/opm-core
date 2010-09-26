G = computeGeometry(cartGrid([30,30,1]));
src = [];
src = addSource(src, 1, 1);
src = addSource(src, G.cells.num, -1);
bc = [];
rock.perm = ones(G.cells.num, 1);
rock.poro = ones(G.cells.num, 1);

W = addWell([], G, rock, 1, 'type', 'bhp', 'val', 200*barsa);

x = initResSol(G, 0, 0);
x = mex_ifsh(x, G, rock, W, bc, src)

plotCellData(G, x.pressure, 'edgec', 'k', 'edgea', .1, 'facea', .625)
view(3), axis tight, grid on, colorbar southoutside
