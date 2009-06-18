function g = makeTest(dims, drop)
physDims   = [1, 1, 1];
g.cartDims = reshape(dims, 1, []);

%% Make pillars

[X,Y,Z] = ndgrid(linspace(0, physDims(1), dims(1) + 1), ...
                 linspace(0, physDims(2), dims(2) + 1), ...
                 linspace(0, physDims(3), 2));
n = prod(dims(1:2) + 1);

lines = zeros([n, 6]);
lines(:, [1, 4]) = reshape(X(:,:,[1, end]), [n, 2]);
lines(:, [2, 5]) = reshape(Y(:,:,[1, end]), [n, 2]);
lines(:, [3, 6]) = reshape(Z(:,:,[1, end]), [n, 2]);
g.COORD = reshape(lines.', [], 1);


layer_thickness = @(x, y) 0.2+0.4*sin(6*pi*(x+rand)).*sin(2*pi*(y+rand));
bottom_layer    = @(x, y) 0.5+0.02*sin(2*pi*x).*sin(2*pi*y);


[X, Y, Z]  = ndgrid(linspace(0, physDims(1), dims(1) + 1), ...
                    linspace(0, physDims(2), dims(2) + 1), ...
                    zeros(dims(3)+1, 1));

Z(:,:,1) = bottom_layer (X(:,:,1) ./ physDims(1), Y(:,:,1) ./ physDims(2));

for k = 2 : dims(3) + 1,
   xi       = X(:,:,k) ./ physDims(1);
   eta      = Y(:,:,k) ./ physDims(2);
   Z(:,:,k) = Z(:,:,k-1) + max(0, layer_thickness(xi, eta)*physDims(3)/dims(3));
end




%% Assign z-coordinates
%ind = @(d) [1, rldecode(2:dims(d), repmat(2,[1,dims(d)-1]), 2), dims(d)+1];
ind = @(d) [1, kron((2:dims(d)), [1,1]), dims(d)+1];
z   = Z(ind(1), ind(2), ind(3));

%% Add fault
%if nargin > 1, z(end/2+1:end,:,:) = z(end/2+1:end,:,:) + drop; end
if nargin > 1, z(end/2+1:end,[1:2:end,end],:) =bsxfun(@plus, z(end/2+1:end,[1:2:end,end],:) , drop(linspace(0,physDims(2), dims(2)+1)));end
if nargin > 1, z(end/2+1:end,2:2:end-1,:) =z(end/2+1:end,3:2:end,:);end
g.ZCORN = z(:);

%% Assign active cells
actnum = reshape(ones(dims, 'int32'), [], 1);

%[i,j,k] = ndgrid(1:dims(1), 1:dims(2), 1:dims(3));
%actnum((i-dims(1)/2).^2+(j-dims(2)/2).^2> 4*min(dims(1)/2, dims(2)/2).^2)=0;
z = reshape(g.ZCORN, 2*g.cartDims);

%actnum(z(1:2:end, 1:2:end,1:2:end)<0.5)=0;

ix = max(1, round(rand(1,1)*prod(dims)));
actnum(ix)=0;
g.ACTNUM = actnum(:);
