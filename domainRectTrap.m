% X1 = (x0, x1) : x-coordinates of lower left and upper right corner
% X2 = (y0, y1) : y-coordinates of lower left and upper right corner or
%                 2 x (nx+1) matrix of y-coordinates
% numElem [= (nx, ny)] : element count in each direction (possibly uniform)
function g = domainRectTrap(X1, X2, numElem)
%% Check input arguments and expand them, if necessary.
validateattributes(numElem, {'numeric'}, {'nonnegative'}, mfilename, 'numElem')
if length(numElem) == 1
  g.numElem = [numElem, numElem];
else
  validateattributes(numElem, {'numeric'}, {'numel', 2}, mfilename, 'numElem')
  g.numElem = reshape(numElem, 1, 2);
end % if
%
validateattributes(X1, {'numeric'}, {'numel', 2}, mfilename, 'X1')
%
if length(X1) == 2
  validateattributes(X1, {'numeric'}, {'numel', 2}, mfilename, 'X1')
  dX1 = (X1(2) - X1(1)) / g.numElem(1);
  X1 = X1(1) : dX1 : X1(2);
else
  validateattributes(X1, {'numeric'},  {'size', [1, g.numElem(1)+1]}, mfilename, 'X1')
end % if
%
if length(X2) == 2
  validateattributes(X2, {'numeric'}, {'numel', 2}, mfilename, 'X2')
  X2 = kron(ones(1, g.numElem(1) + 1), reshape(X2, 2, 1));
else
  validateattributes(X2, {'numeric'}, {'size', [2, g.numElem(1)+1]}, mfilename, 'X2')
end % if
%% Mesh entity counts
g.numV = (g.numElem(1) + 1) * (g.numElem(2) + 1);
g.numT = g.numElem(1) * g.numElem(2);
g.numE = g.numElem(1) * (g.numElem(2) + 1) + g.numElem(2) * (g.numElem(1) + 1);
%% Vertex coordinates (coordV)
g.coordV = zeros(g.numV, 2);
g.coordV(:,1) = repmat( X1.', g.numElem(2) + 1, 1);
g.coordV(:,2) = kron(ones(1, g.numElem(2) + 1), X2(1,:)).' + ...
                kron(0:g.numElem(2), (X2(2,:) - X2(1,:)) / g.numElem(2)).';
%% Mapping Trapezoid -> Vertex (V0T)
g.V0T = zeros(g.numT, 4);
g.V0T(:,1) = kron(0:g.numElem(2)-1, (g.numElem(1) + 1) * ones(1, g.numElem(1))).' + ...
             repmat(1:g.numElem(1), 1, g.numElem(2)).'; % lower left
g.V0T(:,2) = g.V0T(:,1) + 1; % lower right
g.V0T(:,3) = g.V0T(:,2) + g.numElem(1) + 1; % upper right
g.V0T(:,4) = g.V0T(:,3) - 1; % upper left
%% Mapping Trapezoid -> Edge (E0T)
g.E0T = zeros(g.numT, 4);
g.E0T(:,1) = 1 : g.numT; % lower edge
g.E0T(:,2) = g.E0T(:,1) + g.numElem(1); % upper edge
g.E0T(:,3) = g.E0T(:,1) + g.numElem(1) * (g.numElem(2) + 1) + 1; % right edge
g.E0T(:,4) = g.E0T(:,1) + g.numElem(1) * (g.numElem(2) + 1); % left edge
g.E0T(g.numElem(1) : g.numElem(1) : end, 3) = g.numE - g.numElem(2) + 1 : g.numE; % correct right boundary edges
%% Mapping Edge -> Vertex (V0E)
g.V0E = zeros(g.numE, 2);
g.V0E(1:g.numT, :) = g.V0T(:, 1:2); % lower edges in all elements
g.V0E(g.numT:g.numT + g.numElem(1),:) = g.V0T((g.numElem(2) - 1) *  g.numElem(1) : end, [4 3]); % upper edges at top boundary
g.V0E(g.numElem(1) * (g.numElem(2) + 1) + 1 : g.numE - g.numElem(2), :) = g.V0T(:,[1 4]); % left edges
g.V0E(g.numE - g.numElem(2) + 1 : end, :) = g.V0T(g.numElem(1) : g.numElem(1) : end, [2 3]); % right edges
%% Edge lengths and normals (areaE, areaE0T, nuE)
vecE = g.coordV(g.V0E(:, 2), :) - g.coordV(g.V0E(:, 1), :);
areaE = sqrt(vecE(:, 1).^2 + vecE(:, 2).^2);
nuE = vecE * [0, -1; 1, 0] ./ areaE(:, [1, 1]);
g.areaE0T = areaE(g.E0T);
g.nuE0T = zeros(g.numT, 4, 2);
g.nuE0T(:,1,:) = nuE(g.E0T(:,1), :);
g.nuE0T(:,2,:) = -nuE(g.E0T(:,2), :);
g.nuE0T(:,3,:) = nuE(g.E0T(:,3), :);
g.nuE0T(:,4,:) = -nuE(g.E0T(:,4), :);
%% Mapping of neighbouring edges (markE0TE0T)
g.markE0TE0T = cell(1,4);
mapE0E = [2 1 4 3];
for n = 1 : 4
  g.markE0TE0T{n} = sparse(bsxfun(@eq, g.E0T(:,n), g.E0T(:,mapE0E(n))'));
end % for
%% Edge IDs (idE, idE0T)
g.idE = zeros(g.numE, 1);
g.idE(1 : g.numElem(1)) = 1; % Bottom boundary
g.idE(g.numE - g.numElem(2) + 1 : end) = 2; % Right boundary
g.idE(g.numT + 1 : g.numT + g.numElem(1)) = 3; % Top boundary
g.idE(g.numT + g.numElem(1) + 1 : g.numElem(1) : g.numT + g.numElem(1) * (g.numElem(2) + 1)) = 4; % Left boundary
g.idE0T = g.idE(g.E0T);
%% Element-local vertex coordinates (coordV0T)
g.coordV0T = zeros(g.numT, 4, 2);
for k = 1 : 4
  g.coordV0T(:, k, :) = g.coordV(g.V0T(:, k), :);
end % for
%% Element centroids (baryT)
g.baryT = squeeze(sum(g.coordV0T, 2)) / 4;
%% Edge centroids (baryE, baryE0T)
g.baryE = 0.5 * (g.coordV(g.V0E(:,1),:) + g.coordV(g.V0E(:, 2),:));
g.baryE0T = zeros(g.numT, 4, 2);
for k = 1 : 4
  g.baryE0T(:, k, :) = squeeze(g.baryE(g.E0T(:,k), :));
end % for
%% Jacobian of the mapping (J0T) and its determinant (detJ0T), 
% split into constant (1), x1- (2) and x2-contributions (3)
g.J0T = { zeros(g.numT, 2, 2), zeros(g.numT, 2, 2), zeros(g.numT, 2, 2) };
g.J0T{1}(:,1,1) = g.coordV0T(:,2,1) - g.coordV0T(:,1,1);
g.J0T{1}(:,2,1) = g.coordV0T(:,2,2) - g.coordV0T(:,1,2);
g.J0T{1}(:,2,2) = g.coordV0T(:,4,2) - g.coordV0T(:,1,2);
g.J0T{2}(:,2,2) = (g.coordV0T(:,3,2) - g.coordV0T(:,2,2)) - (g.coordV0T(:,4,2) - g.coordV0T(:,1,2));
g.J0T{3}(:,2,1) = (g.coordV0T(:,3,2) - g.coordV0T(:,2,2)) - (g.coordV0T(:,4,2) - g.coordV0T(:,1,2));
g.detJ0T = cell(1,3);
for s = 1 : 3
  g.detJ0T{s} = g.J0T{s}(:,1,1) .* g.J0T{s}(:,2,2) - g.J0T{s}(:,2,1) .* g.J0T{s}(:,1,2);
end % for
%% Element area (areaT)
g.areaT = 0.5 * (g.areaE0T(:,3) + g.areaE0T(:,4)) .* g.J0T{1}(:,1,1);
%% Mapping from reference element to physical element (mapRef2Phy)
g.mapRef2Phy = @(i,X1,X2) g.J0T{1}(:,i,1) * X1 + g.J0T{1}(:,i,2) * X2 + ...
                          g.J0T{2}(:,i,2) * (X1 .* X2) + g.coordV0T(:,1,i) * ones(size(X1));
%% 1D mesh for surface
g.g1D.numT = g.numElem(1);
g.g1D.numV = g.g1D.numT + 1;
g.g1D.idxE2D0T = g.numT + 1 : g.numT + g.numElem(1);
g.g1D.idxV2D0V = g.numElem(2) * (g.numElem(1) + 1) + 1 : (g.numElem(2) + 1) * (g.numElem(1) + 1);
g.g1D.idxT2D0T = bsxfun(@plus, (1 : g.numElem(1) : g.numT).', 0 : g.numElem(2) - 1 ).';
g.g1D.V0T = [(1 : g.g1D.numT).', (2 : g.g1D.numT+1).'];
g.g1D.coordV = [X1.', X2(2,:).'];
g.g1D.coordV0T = zeros(g.g1D.numT, 2, 2);
for k = 1 : 2
  g.g1D.coordV0T(:, k, :) = g.g1D.coordV(g.g1D.V0T(:, k), :);
end % for
g.g1D.idV = zeros(g.g1D.numV, 1); g.g1D.idV(1) = 4; g.g1D.idV(end) = 2;
g.g1D.idV0T = g.g1D.idV(g.g1D.V0T);
g.g1D.nuV0T = repmat([-1 1], g.g1D.numT, 1);
g.g1D.areaT = repmat(dX1, g.g1D.numT, 1);
g.g1D.detJ0T = g.g1D.areaT;
g.g1D.mapRef2Phy = @(X) g.g1D.detJ0T * X + g.g1D.coordV0T(:,1) * ones(size(X));
[c,~,r] = find(g.g1D.idxT2D0T);
g.g1D.markT2DT = sparse(r,c,true(size(r)));
%% Mapping of neighbouring elements (markV0TV0T, markV0TE0T)
g.g1D.markV0TV0T = cell(1,2);
g.g1D.markV0TE0T = cell(1,2);
for n = 1 : 2
  g.g1D.markV0TV0T{n} = sparse(bsxfun(@eq, g.g1D.V0T(:,n), g.g1D.V0T(:,3-n)'));
  g.g1D.markV0TE0T{n} = (double(g.g1D.markT2DT.') * g.markE0TE0T{5-n}) > 0;
end % for
end % function

