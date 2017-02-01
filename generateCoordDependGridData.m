function g = generateCoordDependGridData(g)
%% Edge lengths and normals (areaE, areaE0T, nuE)
vecE = g.coordV(g.V0E(:, 2), :) - g.coordV(g.V0E(:, 1), :);
areaE = sqrt(vecE(:, 1).^2 + vecE(:, 2).^2);
nuE = vecE * [0, -1; 1, 0] ./ areaE(:, [1, 1]);
g.areaE0T = areaE(g.E0T);
g.nuE0T = zeros(g.numT, 4, 2);
g.nuE0T(:, 1, :) = nuE(g.E0T(:, 1), :);
g.nuE0T(:, 2, :) = -nuE(g.E0T(:, 2), :);
g.nuE0T(:, 3, :) = nuE(g.E0T(:, 3), :);
g.nuE0T(:, 4, :) = -nuE(g.E0T(:, 4), :);
%% Element centroids (baryT)
g.baryT = squeeze(sum(g.coordV0T, 2)) / 4;
%% Edge centroids (baryE, baryE0T)
g.baryE = 0.5 * (g.coordV(g.V0E(:, 1), :) + g.coordV(g.V0E(:, 2), :));
g.baryE0T = zeros(g.numT, 4, 2);
for n = 1 : 4
  g.baryE0T(:, n, :) = squeeze(g.baryE(g.E0T(:, n), :));
end % for
%% Jacobian of the mapping (J0T) and its determinant (detJ0T), 
% split into constant (1), x1- (2) and x2-contributions (3)
g.J0T = { zeros(g.numT, 2, 2), zeros(g.numT, 2, 2), zeros(g.numT, 2, 2) };
g.J0T{1}(:, 1, 1) = g.coordV0T(:, 2, 1) - g.coordV0T(:, 1, 1);
g.J0T{1}(:, 2, 1) = g.coordV0T(:, 2, 2) - g.coordV0T(:, 1, 2);
g.J0T{1}(:, 2, 2) = g.coordV0T(:, 4, 2) - g.coordV0T(:, 1, 2);
g.J0T{2}(:, 2, 2) = (g.coordV0T(:, 3, 2) - g.coordV0T(:, 2, 2)) - ...
                    (g.coordV0T(:, 4, 2) - g.coordV0T(:, 1, 2));
g.J0T{3}(:, 2, 1) = (g.coordV0T(:, 3, 2) - g.coordV0T(:, 2, 2)) - ...
                    (g.coordV0T(:, 4, 2) - g.coordV0T(:, 1, 2));
g.detJ0T = cell(1,3);
for s = 1 : 3
  g.detJ0T{s} = g.J0T{s}(:,1,1) .* g.J0T{s}(:,2,2) - g.J0T{s}(:,2,1) .* g.J0T{s}(:,1,2);
end % for
%% Element area (areaT)
g.areaT = 0.5 * (g.areaE0T(:, 3) + g.areaE0T(:, 4)) .* g.J0T{1}(:, 1, 1);
%% Mapping from reference element to physical element (mapRef2Phy)
g.mapRef2Phy = @(i,X1,X2) g.J0T{1}(:,i,1) * X1 + g.J0T{1}(:,i,2) * X2 + ...
                          g.J0T{2}(:,i,2) * (X1 .* X2) + g.coordV0T(:,1,i) * ones(size(X1));
end

