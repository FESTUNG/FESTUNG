function g1D = generateGridData1D(X1, X2, numElem, g2D)
%% Check input arguments and expand them, if necessary.
validateattributes(numElem, {'numeric'}, {'nonnegative', 'numel', 1}, mfilename, 'numElem')
%
validateattributes(X1, {'numeric'}, {'numel', 2}, mfilename, 'X1')
%
if length(X1) == 2
  validateattributes(X1, {'numeric'}, {'numel', 2}, mfilename, 'X1')
  dX1 = (X1(2) - X1(1)) / numElem;
  X1 = X1(1) : dX1 : X1(2);
else
  validateattributes(X1, {'numeric'}, {'numel', numElem+1}, mfilename, 'X1')
end % if
%
if length(X2) == 1
  validateattributes(X2, {'numeric'}, {'numel', 1}, mfilename, 'X2')
  X2 = repmat(X2, 1, numElem+1);
else
  validateattributes(X2, {'numeric'}, {'numel', numElem+1}, mfilename, 'X2')
end % if
%% Mesh entity counts
g1D.numT = numElem;
g1D.numV = g1D.numT + 1;
%% Mapping element -> vertex (V0T)
g1D.V0T = [(1 : g1D.numT).', (2 : g1D.numT+1).'];
%% Vertex coordinates
g1D.coordV = [X1(:), X2(:)];
g1D.coordV0T = zeros(g1D.numT, 2, 2);
for k = 1 : 2
  g1D.coordV0T(:, k, :) = g1D.coordV(g1D.V0T(:, k), :);
end % for
%% Vertex IDs
g1D.idV = zeros(g1D.numV, 1); g1D.idV(1) = 4; g1D.idV(end) = 2;
g1D.idV0T = g1D.idV(g1D.V0T);
g1D.nuV0T = repmat([-1 1], g1D.numT, 1);
%% Element sizes
g1D.areaT = g1D.coordV0T(:, 2) - g1D.coordV0T(:, 1);
%% Mapping to reference element and its Jacobian
g1D.detJ0T = g1D.areaT;
g1D.mapRef2Phy = @(X) g1D.detJ0T * X + g1D.coordV0T(:, 1) * ones(size(X));
%% Mapping of neighbouring elements (markV0TV0T)
g1D.markV0TV0T = cell(1,2);
for n = 1 : 2
  g1D.markV0TV0T{n} = sparse(bsxfun(@eq, g1D.V0T(:,n), g1D.V0T(:,3-n)'));
end % for
%% Create connectivity between 2D and 1D mesh, if given.
if nargin == 4
  validateattributes(g2D, {'struct'}, {}, mfilename, 'g2D')
  assert(g1D.numT == g2D.numElem(1), 'Number of horizontal 2D elements does not match given number');
  g1D.idxE2D0T = g2D.numT + 1 : g2D.numT + g2D.numElem(1);
  g1D.idxV2D0V = g2D.numElem(2) * (g2D.numElem(1) + 1) + 1 : (g2D.numElem(2) + 1) * (g2D.numElem(1) + 1);
  g1D.idxT2D0T = bsxfun(@plus, (1 : g2D.numElem(1) : g2D.numT).', 0 : g2D.numElem(1) - 1 ).';
  [c, ~, r] = find(g1D.idxT2D0T);
  g1D.markT2DT = sparse(r, c, true(size(r)));
  g1D.markV0TE0T = cell(1,2);
  for n = 1 : 2
    g1D.markV0TE0T{n} = (double(g1D.markT2DT.') * g2D.markE0TE0T{5-n}) > 0;
  end % for
end % if
end

