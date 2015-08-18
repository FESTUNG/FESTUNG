% Assembles lists describing the geometric and topological
% properties of a triangulation.
%
%===============================================================================
%> @file generateGridData.m
%>
%> @brief Assembles lists describing the geometric and topological
%>        properties of a triangulation.
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
%>   - <code>B</code> @f$[\#\mathcal{T} \times 2 \times 2]@f$:
%>     transformation matrices @f$\mathsf{{B}}_k@f$
%>   - <code>areaE</code> @f$[\#\mathcal{E} \times 1]@f$:
%>     edge lengths @f$E_{n}@f$
%>   - <code>areaEOT</code> @f$[\#\mathcal{T} \times 3]@f$:
%>     edge lengths @f$E_{kn}@f$
%>   - <code>areaT</code> @f$[\#\mathcal{T} \times 1]@f$:
%>     triangle areas @f$|\mathcal{T}|@f$
%>   - <code>coordV</code> @f$[\#\mathcal{V} \times 2]@f$:
%>     vertex coordinates @f$\mathbf{a}_{n}@f$
%>   - <code>coordV0T</code> @f$[\#\mathcal{T} \times 3 \times 2]@f$:
%>     vertex coordinates @f$\mathbf{a}_{kn}@f$
%>   - <code>nuE</code> @f$[\#\mathcal{E} \times 2]@f$:
%>     global edge normals @f$\mathbf{\nu}_n@f$
%>   - <code>nuE0T</code> @f$[\#\mathcal{T} \times 3 \times 2]@f$:
%>     local edge normals @f$\mathbf{\nu}_{kn}@f$, exterior to @f$T_k@f$
%>   - <code>baryT</code> @f$[\#\mathcal{T} \times 2]@f$:
%>     coordinates of the barycenter of each triangle @f$T_k@f$
%>   - <code>baryE</code> @f$[\#\mathcal{E} \times 2]@f$:
%>     coordinates of the midpoint of each edge @f$E_n@f$
%>   - <code>baryE0T</code> @f$[\#\mathcal{T} \times 3 \times 2]@f$:
%>     local edge midpoints of each edge @f$E_{kn}@f$
%> - Topological data:
%>   - <code>E0T</code> @f$[\#\mathcal{T} \times 3]@f$:
%>     global edge indices of triangles
%>   - <code>idE0T</code> @f$[\#\mathcal{T} \times 3]@f$:
%>     edge IDs for edges @f$E_{kn}@f$ (used to identify the interior and 
%>     boundary edges as well as Dirichlet and Neumann edges)
%>   - <code>markE0TE0T</code> @f$[3 \times 3 \text{ (cell)}]@f$:
%>     the @f$(n^-, n^+)@f$th entry of this cell is a sparse @f$K \times K@f$
%>     array whose @f$(k^-,k^+)@f$th entry is one if @f$E_{k^-n^-}=E_{k^+n^+}@f$
%>   - <code>T0E</code> @f$[\#\mathcal{E} \times 2]@f$:
%>     global indices of the triangles sharing edge @f$E_n@f$ in the order
%>     dictated by the direction of the global normal on @f$E_n@f$
%>     (i.e. @f$T^-, T^+ \text{ if } \mathbf{\nu}_{E_n}=\mathbf{\nu}_{T^-}@f$)
%>   - <code>V0E</code> @f$[\#\mathcal{E} \times 2]@f$:
%>     global indices of the vertices sharing edge @f$E_n@f$ ordered according
%>     to the global edge orientation (the latter is given by rotating
%>     counter-clockwise by @f$\pi/2@f$ the global edge normal to @f$E_n@f$)
%>   - <code>V0T</code> @f$[\#\mathcal{T} \times 3]@f$:
%>     global vertex indices of triangles accounting for the counter-clockwise
%>     ordering
%>   - <code>V2T</code> @f$[\#\mathcal{V} \times 3]@f$:
%>     global vertex indices of triangles accounting for the counter-clockwise
%>     ordering
%>   - <code>V2E</code> @f$[???]@f$:
%>     ???
%>   - <code>E0E</code> @f$[\#\mathcal{E} \times 2]@f$:
%>     ???
%>
%> Furthermore, it provides the following counts of mesh entities:
%> 
%> - <code>numT</code> (scalar): number of triangles @f$\#\mathcal{T}=K@f$
%> - <code>numE</code> (scalar): number of edges @f$\#\mathcal{E}@f$
%> - <code>numV</code> (scalar): number of vertices @f$\#\mathcal{V}@f$
%>
%> Note that the lists <code>g.idE</code> and <code>g.idE0T</code>
%> storing the global and local edge indices are <i>not</i> generated and
%> have to be defined manually after calling <code>generateGridData</code>.
%>
%> @par Example
%> @parblock
%> @code
%> g = generateGridData([0, -1; sqrt(3), 0; 0, 1; -sqrt(3), 0], [4,1,3; 1,2,3]);
%> g.idE = (abs(g.nuE(:,2)) > 0) .* ((g.nuE(:,1)>0) + (g.nuE(:,2)>0)*2+1);
%> visualizeGrid(g)
%> @endcode
%>
%> produces the following output:
%> @image html  generateGridData.png
%> @endparblock
%>
%> @param  coordV  Contains the @f$x^1@f$ and @f$x^2@f$ coordinates of the 
%>                 grid vertices (using a <i>global</i> index) 
%>                 @f$[\#\mathcal{V}\times 2]@f$
%> @param  V0T     The global vertex indices for each triangle with a 
%>                 counter-clockwise ordering @f$[\#\mathcal{T}_h \times 3]@f$
%> @retval g       A struct containing the lists that describe the triangulation
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2015 Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function g = generateGridData(coordV, V0T)
g.coordV = coordV;
g.V0T = V0T;
g.numT = size(g.V0T, 1);
g.numV = size(g.coordV, 1);
% The following implicitely defines the signs of the edges.
g.V2T  = sparse(g.V0T(:, [1 2 3 1 2 3 1 2 3]), g.V0T(:, [2 3 1 2 3 1 2 3 1]), ...
  [(1:g.numT)',zeros(g.numT,3),(1:g.numT)',zeros(g.numT,3),(1:g.numT)'],g.numV,g.numV);
% The following implicitely defines the edge numbers.
[r, c] = find(triu(g.V2T + g.V2T'));
g.V2E = sparse(r, c, 1 : size(r, 1), g.numV, g.numV);
g.V2E = g.V2E + g.V2E';
idxE = full(g.V2E(sub2ind([g.numV,g.numV],g.V0T(end:-1:1,[1,2,3]),g.V0T(end:-1:1,[2,3,1]))))';
g.V0E(idxE(:), 1) = reshape(g.V0T(end:-1:1, [1,2,3])', 3*g.numT, 1);
g.V0E(idxE(:), 2) = reshape(g.V0T(end:-1:1, [2,3,1])', 3*g.numT, 1);
g.T0E(idxE(:), 1) = reshape(full(g.V2T(sub2ind([g.numV,g.numV], ...
  g.V0T(end:-1:1,[1,2,3]), g.V0T(end:-1:1,[2,3,1]))))', 3*g.numT, 1);
g.T0E(idxE(:), 2) = reshape(full(g.V2T(sub2ind([g.numV,g.numV], ...
  g.V0T(end:-1:1,[2,3,1]), g.V0T(end:-1:1,[1,2,3]))))', 3*g.numT, 1);
g.numE = size(g.V0E, 1);
vecE = g.coordV(g.V0E(:, 2), :) - g.coordV(g.V0E(:, 1), :);
g.areaE = (vecE(:, 1).^2 + vecE(:, 2).^2).^(1/2);
g.nuE = vecE * [0,-1; 1,0] ./ g.areaE(:, [1, 1]);
g.areaT = ...
  ( g.coordV(g.V0T(:,1),1).*g.coordV(g.V0T(:,2),2) + g.coordV(g.V0T(:,2),1).*g.coordV(g.V0T(:,3),2) ...
  + g.coordV(g.V0T(:,3),1).*g.coordV(g.V0T(:,1),2) - g.coordV(g.V0T(:,1),1).*g.coordV(g.V0T(:,3),2) ...
  - g.coordV(g.V0T(:,2),1).*g.coordV(g.V0T(:,1),2) - g.coordV(g.V0T(:,3),1).*g.coordV(g.V0T(:,2),2) )/2;
g.baryT = (g.coordV(g.V0T(:,1),:)+g.coordV(g.V0T(:,2),:)+g.coordV(g.V0T(:,3),:))/3;
g.E0T = full(g.V2E(sub2ind([g.numV,g.numV],g.V0T(:,[2,3,1]),g.V0T(:,[3,1,2]))));
g.areaE0T = g.areaE(g.E0T);
sigE0T = 1-2*(bsxfun(@eq, reshape(g.T0E(g.E0T,2),g.numT,3),(1:g.numT)'));
g.baryE = 0.5 * (g.coordV(g.V0E(:, 1), :) + g.coordV(g.V0E(:, 2), :));
for n = 1 : 3
  for m = 1 : 2
    g.coordV0T(:, n, m) = g.coordV(g.V0T(:, n), m)';
    g.baryE0T(:, n, m) = g.baryE(g.E0T(:, n), m)';
    g.nuE0T(:, n, m) = g.nuE(g.E0T(:, n), m)'.* sigE0T(:, n)';
  end 
  Tn = sigE0T(:, n) == 1;  Tp = ~Tn;
  g.E0E(g.E0T(Tn, n), 1) = n;  g.E0E(g.E0T(Tp, n), 2) = n;
end % for
for m = 1 : 2
  g.B(:, m, 1) = g.coordV0T(:, 2, m) - g.coordV0T(:, 1, m);
  g.B(:, m, 2) = g.coordV0T(:, 3, m) - g.coordV0T(:, 1, m);
end % for
markEint = g.E0E(:, 2) ~= 0; % mark interior edges
g.markE0TE0T = cell(3, 3);
for nn = 1 : 3
  for np = 1 : 3
    g.markE0TE0T{nn,np} = sparse(g.numT, g.numT);
    markEn = g.E0E(:, 1) == nn;  markEp = g.E0E(:, 2) == np;
    idx = markEn & markEp & markEint;
    g.markE0TE0T{nn, np}(sub2ind([g.numT, g.numT], g.T0E(idx, 1), g.T0E(idx, 2))) = 1;
    markEn = g.E0E(:, 2) == nn;  markEp = g.E0E(:, 1) == np;
    idx = markEn & markEp & markEint;
    g.markE0TE0T{nn, np}(sub2ind([g.numT, g.numT], g.T0E(idx, 2), g.T0E(idx, 1))) = 1;
  end % for
end % for
end % function
