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
g.V0E(g.numT + 1:g.numT + g.numElem(1),:) = g.V0T((g.numElem(2) - 1) *  g.numElem(1) + 1: end, [4 3]); % upper edges at top boundary
g.V0E(g.numElem(1) * (g.numElem(2) + 1) + 1 : g.numE - g.numElem(2), :) = g.V0T(:,[1 4]); % left edges
g.V0E(g.numE - g.numElem(2) + 1 : end, :) = g.V0T(g.numElem(1) : g.numElem(1) : end, [2 3]); % right edges
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
%% Generate grid data that changes when coordV, coordV0T changes
g.generateCoordDependGridData = getFunctionHandle('darcyVert/generateCoordDependGridData');
g = g.generateCoordDependGridData(g);
end % function

