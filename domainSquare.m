% Generates a mesh for the unit square.
%
%===============================================================================
%> @file domainSquare.m
%>
%> @brief Generates a mesh for the unit square.
%===============================================================================
%>
%> @brief Generates a mesh for the unit square.
%>
%> This function allows to generate a Friedrich-Keller triangulation for 
%> the unit square @f$\{(0,0), (1,0), (1,1), (0,1)\}@f$
%>
%> @par Example
%> @parblock
%> @code
%> g = domainSquare(1/3);
%> visualizeGrid(g)
%> @endcode
%> produces the following output:
%> @image html  domainSquare.png
%> @endparblock
%>
%> @param  h      The maximum diameter of an element.
%> @retval g      A struct containing the lists that describe the
%>                triangulation, as explained in <code>generateGridData()</code>.
%>                <code>idE</code> and <code>ideE0T</code> are filled correctly.
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
function g = domainSquare(h)
dim = ceil(1/h); % number of edges per side of the unit square
h = 1/dim;
%% Build coordV.
[X, Y] = meshgrid(0:h:1);
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
g.idE(g.baryE(:, 2) == 0) = 1; % south
g.idE(g.baryE(:, 1) == 1) = 2; % east
g.idE(g.baryE(:, 2) == 1) = 3; % north
g.idE(g.baryE(:, 1) == 0) = 4; % west
g.idE0T = g.idE(g.E0T); % local edge IDs
end % function
