% Generates a mesh for a polygonal domain with given mesh width.

%===============================================================================
%> @file
%>
%> @brief Generates a mesh for a polygonal domain with given mesh width.
%===============================================================================
%>
%> @brief Generates a mesh for a polygonal domain with given mesh width.
%>
%> This function allows to generate a triangulation for arbitrary polygonal
%> domains.
%>
%> @par Example
%> @parblock
%> @code
%> X1 = [0.0, 0.5, 0.5, 1.0, 1.0, 0.0];
%> X2 = [0.0, 0.0, 0.5, 0.5, 1.0, 1.0];
%> g = domainPolygon(X1, X2, 0.5);
%> visualizeGrid(g)
%> @endcode
%> produces the following output:
%> @image html  domainPolygon.png
%> @endparblock
%>
%> @param  X1,X2  The @f$x^1,x^2@f$ coordinates of the boundary points 
%>                that describe the domain.
%> @param  h      The maximum diameter of an element.
%> @retval g      A struct containing the lists that describe the
%>                triangulation, as explained in <code>generateGridData()</code>.
%>                <code>idE</code> and <code>ideE0T</code> are filled correctly.
%>
%> @note  MATLABs Partial Differential Equation Toolbox or at least the
%>        function <code>initMesh()</code> is required for this routine.
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
function g = domainPolygon(X1, X2, h)
assert(length(X1) >= 3, 'At least 3 points are required for a 2D domain')
assert(isequal(size(X1), size(X2)), 'X1 and X2 must be of same size')
gd = [2; length(X1(:)); X1(:); X2(:)]; % geometry description
sf = 'polygon';                        % set formula
ns = double('polygon')';               % name space
[p, e, t] = initmesh(decsg(gd,sf,ns), 'Hmax', h);
g  = generateGridData(p', t(1:3, :)');
g.idE = zeros(g.numE, 1);
g.idE(g.V2E(sub2ind(size(g.V2E),e(1,:),e(2,:)))) = e(5,:);
g.idE0T = g.idE(g.E0T); % local edge IDs
end % function
