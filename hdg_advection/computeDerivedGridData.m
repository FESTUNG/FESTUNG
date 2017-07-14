% Computes various additional fields needed throughout the advection
% problem presented in @ref RAWFK2016 . These fields only contain
% information based on the underlying grid and are therefore constant in
% time but used in every step. They are saved as part of the grid.

%===============================================================================
%> @file hdg_advection/computeDerivedGridData.m
%>
%> @brief Computes various additional fields needed throughout the advection
%>        problem presented in @ref RAWFK2016 . These fields only contain
%>        information based on the underlying grid and are therefore constant in
%>        time but used in every step. They are saved as part of the grid.
%===============================================================================
%>
%> @brief Computes various additional fields needed throughout the advection
%>        problem presented in @ref RAWFK2016 . These fields only contain
%>        information based on the underlying grid and are therefore constant in
%>        time but used in every step. They are saved as part of the grid.
%>
%> For every vertex of a triangle the columns of g.markV0TT0V indicate which triangle 
%> corresponds to this vertex.
%>
%> g.areaE0TbdrD provides the elementwise products of the vectors 
%> <code>g.areaE0T</code>, and <code>markE0TbdrD</code>.
%> @f$[3 \text{ cell}]@f$
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @retval g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2015 Hennes Hajduk, Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function g = computeDerivedGridData(g)
g.markSideE0T = true(g.numT, 3, 2);
for n = 1 : 3
  % Mark element to have local id 1 or 2 at the given edge
  g.markSideE0T(:, n, 1) = g.T0E(g.E0T(:, n), 2) ~= (1:g.numT)';
  g.markSideE0T(:, n, 2) = ~g.markSideE0T(:, n, 1);
end % for
end % function