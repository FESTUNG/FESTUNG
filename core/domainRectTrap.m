% Generates a trapezoidal mesh for a polygonal domain with given number of
% elements per dimension.

%===============================================================================
%> @file
%>
%> @brief Generates a trapezoidal mesh for a polygonal domain with given 
%>        number of elements per dimension.
%===============================================================================
%>
%> @brief Generates a trapezoidal mesh for a polygonal domain with given 
%>        number of elements per dimension.
%>
%> This function allows to generate a trapezoidal mesh for (almost) arbitrary 
%> polygonal domains. The domain shape is restricted by the fact that
%> lateral edges are always vertical, i.e., aligned with the x2-axis (thus,
%> producing trapezoidal elements).
%>
%> The shape of the top/bottom boundary can be specified as point-wise 
%> coordinates. For rectangular domains it is sufficient to specify the 
%> coordinates of the lower-left and upper-right corners.
%>
%> The number of elements @f$N@f$ per dimension is specified as a function 
%> argument, thus the resulting mesh has @f$N^2@f$ elements. Optionally,
%> different element counts can be specified for x1- and x2-direction.
%>
%> The resulting mesh can be visualized using the function
%> <code>visualizeGridQuadri()</code>.
%>
%> The bottom boundary is assigned id 1, right boundary has id 2, top
%> boundary has id 3, and left boundary is assigned id 4.
%>
%> @par Example
%> @parblock
%> @code
%> X1 = [0.0, 100.0];
%> X2 = [0.0, 1.0];
%> g = domainRectTrap(X1, X2, 10);
%> @endcode
%> Creates a mesh for the domain @f$[0, 100] \times [0,1]@f$ with 10
%> elements in each direction.
%> @endparblock
%>
%> @par Example
%> @parblock
%> @code
%> X1 = [0.0, 1.0, 2.0, 2.5, 3.0];
%> X2 = [0.0, 0.1, 0.2, 0.1, 0.2; ...
%>       2.0, 1.8, 1.7, 1.6, 1.5];
%> g = domainRectTrap(X1, X2, [4, 3]);
%> visualizeGridQuadri(g);
%> @endcode
%> Creates a mesh for a domain with x1-values from 0 to 3, 
%> bottom-boundary between 0 and 0.2, and top-boundary between 1.5 and 2.
%> The resulting mesh has 12 elements (4 in x1-direction and 3 in
%> x2-direction. Note that the number of coordinates must match the number
%> of vertices in x1-direction (i.e., number of elements plus one).
%> The resulting grid looks like this:
%> @image html  domainRectTrap.png 
%> @endparblock
%>
%> @param  X1       The @f$x^1@f$ coordinates of the boundary points 
%>                  that describe the domain.
%>                  @f$[2 \times 1]@f$ or @f$[N_1 \times 1]@f$
%> @param  X2       The @f$x^2@f$ coordinates of the boundary points 
%>                  that describe the domain.
%>                  @f$[2 \times 1]@f$ or @f$[N_1 \times 2]@f$
%> @param  numElem  The number of elements, either as two-dimensional
%>                  vector @f$[N_1, N_2]@f$ or as scalar @f$N@f$.
%>                  @f$[2 \times 1]@f$ or @f$[1 \times 1]@f$
%> @retval g        A struct containing the lists that describe the
%>                  triangulation.
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
%> 
%> @par License
%> @parblock
%> This program is free software: you can redistribute it and/or modify
%> it under the terms of the GNU General Public License as published by
%> the Free Software Foundation, either version 3 of the License, or
%> (at your option) any later version.
%>
%> This program is distributed in the hope that it will be useful,
%> but WITHOUT ANY WARRANTY; without even the implied warranty of
%> MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%> GNU General Public License for more details.
%>
%> You should have received a copy of the GNU General Public License
%> along with this program.  If not, see <http://www.gnu.org/licenses/>.
%> @endparblock
%
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
for n = 1 : 4
  g.markE0TE0T{n} = sparse(bsxfun(@eq, g.E0T(:, n), g.E0T(:, mapLocalEdgeIndexQuadri(n))'));
end % for
% Fix a bug in GNU Octave 4.0.0's implementation of sparse matrix concatenation
if isOctave
  for n = 1 : 4
    g.markE0TE0T{n} = g.markE0TE0T{n} + 0 * speye(size(g.markE0TE0T{n}, 1), size(g.markE0TE0T{n}, 2));
  end % for n
end % if
%% Mapping of neighbouring vertices (markV0TV0T)
g.markV0TV0T = cell(4, 4);
try
  for nn = 1 : 4
    for np = 1 : 4
      g.markV0TV0T{nn, np} = sparse(bsxfun(@eq, g.V0T(:,nn), g.V0T(:,np)'));
    end % for np
  end % for nn
catch
  g = rmfield(markV0TV0T);
end
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
g.generateCoordDependGridData = @generateCoordDependGridData;
g = g.generateCoordDependGridData(g);
end % function
%
%===============================================================================
%>
%> @brief Helper routine that allows to recreate all coordinate-dependend
%>        data-structures.
%>
%> It can be evaluated using the <code>function_handle</code> in the grid 
%> data structure: <code>g = g.generateCoordDependGridData(g)</code>.
%> Useful when adapting the mesh.
%
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
g.detJ0T = { g.J0T{1}(:, 1, 1) .* g.J0T{1}(:, 2, 2), ...
             g.J0T{1}(:, 1, 1) .* g.J0T{2}(:, 2, 2), ...
             zeros(g.numT, 1) };
%% Element area (areaT)
g.areaT = 0.5 * (g.areaE0T(:, 3) + g.areaE0T(:, 4)) .* g.J0T{1}(:, 1, 1);
%% Mapping from reference element to physical element (mapRef2Phy)
J0T1 = g.J0T{1}; J0T2 = g.J0T{2}; a1 = g.coordV0T(:,1,:);
g.mapRef2Phy = @(i,X1,X2) J0T1(:,i,1) * X1 + J0T1(:,i,2) * X2 + J0T2(:,i,2) * (X1 .* X2) + a1(:,1,i) * ones(size(X1));
end
