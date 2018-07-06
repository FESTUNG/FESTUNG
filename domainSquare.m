% Generates a mesh for a square.

%===============================================================================
%> @file
%>
%> @brief Generates a Friedrichs-Keller triangulation on a square.
%===============================================================================
%>
%> @brief Generates a Friedrichs-Keller triangulation on a square.
%>
%> This function generates a Friedrichs-Keller triangulation on a
%> square-bounded domain @f$\{(a,a), (b,a), (b,b), (a,b)\}@f$.
%> If no bounds for @f$x, y@f$ are specified, it defaults to the unit
%> square @f$[0,1]^2@f$.
%> Optionally, the mesh can be transformed such that every element becomes
%> an equilateral triangle with the same maximum edge length h.  This mesh
%> is a Voronoi grid.  The square-shaped boundary becomes a parallelogram
%> in this case.
%> 
%>
%> @par Example
%> @parblock
%> @code
%> g1 = domainSquare(1/3);
%> visualizeGrid(g1)
%> g2 = domainSquare(1/3, 0, 1, true);
%> visualizeGrid(g2)
%> @endcode
%> produces the following output:
%> @image html  domainSquare.png  Grid g1
%> @image html  domainSquareEquilateral.png   Grid g2 (isEquilateral = true)
%> @endparblock
%>
%> @param  h      The maximum diameter of an element.
%> @param  a     (optional) Lower bound @f$a@f$. Defaults to 0.
%> @param  b     (optional) Upper bound @f$b@f$. Defaults to 1.
%> @param  isEquilateral (optional) Transform to equilateral mesh.
%> @retval g      A struct containing the lists that describe the
%>                triangulation, as explained in <code>generateGridData()</code>.
%>                <code>idE</code> and <code>ideE0T</code> are filled correctly.
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2018 Florian Frank, Balthasar Reuter, Vadym Aizinger
%> @author Balthasar Reuter, 2017.
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
function g = domainSquare(h, a, b, isEquilateral)
if nargin < 2, a = 0; end
if nargin < 3, b = 1; end
dim = ceil((b-a)/h); % number of edges per side of the unit square
h = (b-a)/dim;
%% Build coordV.
[X, Y] = meshgrid(a:h:b);
if nargin == 4 && isEquilateral == true % Transformation to become equilateral.
  xTrans = repmat((0:h/2:(b-a)/2).', 1, dim+1);
  X = X + xTrans;
  yScale = sqrt(3)/2;
  Y = (Y-(b-a))*yScale + (b-a);
end
Xlist = reshape(X, length(X)^2, 1);  Ylist = reshape(Y, length(X)^2, 1);
coordV = [Xlist, Ylist];
%% Build V0T.
pat1 = [1,dim+2,2]; % pattern of "lower-left" triangles
V0T1 = repmat(pat1, dim*(dim+1), 1) + repmat((0:dim*(dim+1)-1)', 1, 3);
V0T1(dim+1 : dim+1 : dim*(dim+1), :) = [];
pat2 = [dim+2,dim+3,2];
V0T2 = repmat(pat2, dim*(dim+1), 1) + repmat((0:dim*(dim+1)-1)', 1, 3);
V0T2(dim+1 : dim+1 : dim*(dim+1), :) = [];
%% Generate grid data and boundary IDs
g = generateGridData(coordV, [V0T1; V0T2]);
g.idE = zeros(g.numE, 1);
g.idE(g.baryE(:, 2) == a) = 1; % south
g.idE(g.baryE(:, 1) == b) = 2; % east
g.idE(g.baryE(:, 2) == b) = 3; % north
g.idE(g.baryE(:, 1) == a) = 4; % west
g.idE0T = g.idE(g.E0T); % local edge IDs
end % function
