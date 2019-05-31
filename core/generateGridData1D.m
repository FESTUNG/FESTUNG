% Assembles lists describing the geometric and topological
% properties of a one-dimensional triangulation.

%===============================================================================
%> @file
%>
%> @brief Assembles lists describing the geometric and topological
%>        properties of a one-dimensional triangulation.
%===============================================================================
%>
%> @brief Assembles lists describing the geometric and topological properties  
%>        of a triangulation @f$\mathcal{T}_h@f$ and stores them in the output 
%>        variable <code>g</code> of type <code>struct</code>.
%>
%> All lists fall in one of two categories: "geometric data" or 
%> "topological data". The provided lists are (dimension in
%> brackets):
%>
%> - Geometric data:
%>   - <code>areaT</code> @f$[\#\mathcal{T} \times 1]@f$:
%>     element sizes @f$|\mathcal{T}|@f$
%>   - <code>coordV</code> @f$[\#\mathcal{V} \times 2]@f$:
%>     vertex coordinates @f$\mathbf{a}_{n}@f$
%>   - <code>coordV0T</code> @f$[\#\mathcal{T} \times 2 \times 2]@f$:
%>     element-local vertex coordinates @f$\mathbf{a}_{kn}@f$
%>   - <code>nuV0T</code> @f$[\#\mathcal{T} \times 2]@f$:
%>     local edge normals @f$\nu_{kn}@f$ (i.e., @f$\pm 1@f$), exterior to @f$T_k@f$
%>   - <code>mapRef2Phy</code> @f$[1 \times 1 \mathtt{function\_handle}]@f$:
%>     maps a given coordinate in the reference interval [0,1] to the
%>     corresponding physical coordinate in all elements.
%>   - <code>detJ0T</code> @f$[\#\mathcal{T} \times 1]@f$:
%>     Determinant of the mappings Jacobian.
%> - Topological data:
%>   - <code>V0T</code> @f$[\#\mathcal{T} \times 2]@f$:
%>     global vertex indices of nodes in elements
%>   - <code>idV</code> @f$[\#\mathcal{V} \times 1]@f$:
%>     IDs for vertices @f$\mathbf{a}_{n}@f$ (used to identify the interior and 
%>     boundary vertices as well as Dirichlet and Neumann vertices)
%>   - <code>idV0T</code> @f$[\#\mathcal{T} \times 2]@f$:
%>     IDs for element-local vertices @f$\mathbf{a}_{kn}@f$ (used to 
%>     identify the interior and boundary vertices)
%>   - <code>markV0TV0T</code> @f$[2 \times 1 \text{ (cell)}]@f$:
%>     the @f$n^-@f$th entry of this cell is a sparse @f$K \times K@f$
%>     array whose @f$(k^-,k^+)@f$th entry is one if 
%>     @f$\mathbf{a}_{k^-n^-}=\mathbf{a}_{k^+n^+}@f$
%>
%> Furthermore, it provides the following counts of mesh entities:
%> 
%> - <code>numT</code> (scalar): number of elements @f$\#\mathcal{T}=K@f$
%> - <code>numV</code> (scalar): number of vertices @f$\#\mathcal{V}@f$. It
%>   holds @f$\#\mathcal{V} = \#\mathcal{T} + 1@f$.
%>
%> @param  X1      Contains the @f$x^1@f$ coordinates of the 
%>                 grid vertices (using a <i>global</i> index) 
%>                 @f$[\#\mathcal{V}\times 1]@f$ or the minimum- and
%>                 maximum @f$x^1@f$ coordinate @f$[2\times 1]@f$
%> @param  X2      Contains the @f$x^2@f$ coordinates of the 
%>                 grid vertices (using a <i>global</i> index) 
%>                 @f$[\#\mathcal{V}\times 1]@f$ or a single @f$x^2@f$ 
%>                 coordinate that is the same for all vertices 
%>                 @f$[1\times 1]@f$
%> @param  numElem The number of elements @f$\#\mathcal{T}@f$ @f$[1\times 1]@f$
%> @param  g2D     (optional) A quadrilateral mesh (e.g., trapezoidal mesh
%>                 as provided by <code>domainRectTrap()</code>) with which
%>                 the one-dimensional mesh is aligned.
%> @retval g       A struct containing the lists that describe the triangulation
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
function g1D = generateGridData1D(X1, X2, numElem, g2D)
%% Check input arguments and expand them, if necessary.
validateattributes(numElem, {'numeric'}, {'nonnegative', 'numel', 1}, mfilename, 'numElem')
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
detJ0T = g1D.detJ0T; a1 = g1D.coordV0T(:, 1);
g1D.mapRef2Phy = @(X) detJ0T * X + a1 * ones(size(X));
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
  % Fix a bug in GNU Octave 4.0.0's implementation of sparse matrix concatenation
  if isOctave
    g1D.markT2DT = g1D.markT2DT + 0 * speye(size(g1D.markT2DT, 1), size(g1D.markT2DT, 2));
    for n = 1 : 2
      g1D.markV0TE0T{n} = g1D.markV0TE0T{n} + 0 * speye(size(g1D.markV0TE0T{n}, 1), size(g1D.markV0TE0T{n}, 2));
    end % for
  end % if
end % if
end

