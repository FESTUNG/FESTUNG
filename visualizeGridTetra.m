% Visualize the quadrilateral grid with global and local indices.

%===============================================================================
%> @file visualizeGridTetra.m
%>
%> @brief Visualize the quadrilateral grid with global and local indices.
%===============================================================================
%>
%> @brief Visualize the quadrilateral grid with global and local indices.
%>
%> @par Example
%> @parblock
%> @code
%> g = domainRectTrap([0, 1], [0, 1], [4, 2]);
%> visualizeGridTetra(g)
%> @endcode
%> produces the following output:
%> @image html  visualizeGridTetra.png
%> @endparblock
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a quadrilateral grid (see, e.g., 
%>                    <code>domainRectTrap()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
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
function visualizeGridTetra(g)
figure('Color', [1, 1, 1]); % white background
hold('on'),  axis('off')
% daspect([1, 1, 1]) % adjust aspect ration, requires Octave >= 3.8
textarray = @(x1,x2,s,col) arrayfun(@(a,b,c) text(a,b,int2str(c),'Color',col), x1, x2, s);
%% Trapezoidal boundaries.
X1 = reshape(g.coordV(:,1), [], g.numElem(2) + 1).';
X2 = reshape(g.coordV(:,2), [], g.numElem(2) + 1).';
surf(X1, X2, zeros(size(X1)), 'facecolor', 'none');
%% Local edge numbers.
w = [4/9, 4/9, 1/18, 1/18; 1/18, 1/18, 4/9, 4/9; ...
     1/18, 4/9, 4/9, 1/18; 4/9, 1/18, 1/18, 4/9].';
for kE = 1 : 4
  textarray(reshape(g.coordV(g.V0T,1),g.numT,4)*w(:,kE), ...
            reshape(g.coordV(g.V0T,2),g.numT,4)*w(:,kE), kE*ones(g.numT, 1), 'blue')
end % for
%% Global vertex numbers.
textarray(g.coordV(:,1), g.coordV(:,2), (1:g.numV)', 'green');
%% Local vertex numbers.
w = ones(4) / 18 + eye(4) * (5/6 - 1/18);
for kE = 1 : 4
  textarray(reshape(g.coordV(g.V0T,1),g.numT,4)*w(:,kE), ...
            reshape(g.coordV(g.V0T,2),g.numT,4)*w(:,kE), kE*ones(g.numT, 1), 'blue')
end % for
%% Global edge numbers.
textarray(g.baryE(:,1), g.baryE(:,2), (1:g.numE)', 'green');
%% Element numbers.
textarray(g.baryT(:,1), g.baryT(:,2), (1:g.numT)', 'green');
%% Edge IDs.
markE0Tbdr = g.idE0T ~= 0;
[r, c] = find(markE0Tbdr);
ind1 = sub2ind([g.numT, 4, 2], r, c, 1 * ones(size(r)));
ind2 = sub2ind([g.numT, 4, 2], r, c, 2 * ones(size(r)));
textarray(g.baryE0T(ind1) + g.nuE0T(ind1).*g.areaE0T(markE0Tbdr)/16, ...
          g.baryE0T(ind2) + g.nuE0T(ind2).*g.areaE0T(markE0Tbdr)/8, ...
          g.idE0T(markE0Tbdr), 'red')
end % function
