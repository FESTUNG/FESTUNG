% Computes an additional field needed throughout the advection
% problem presented in @ref JRASK2017 . This field only contains
% information based on the underlying grid and is therefore constant in
% time but used in every step. It is saved as part of the grid.

%===============================================================================
%> @file hdg_advection/computeDerivedGridData.m
%>
%> @brief Computes an additional field needed throughout the advection
%>        problem presented in @ref JRASK2017 . This field only contains
%>        information based on the underlying grid and is therefore constant in
%>        time but used in every step. It is saved as part of the grid.
%===============================================================================
%>
%> @brief Computes an additional field needed throughout the advection
%>        problem presented in @ref JRASK2017 . This field only contains
%>        information based on the underlying grid and is therefore constant in
%>        time but used in every step. It is saved as part of the grid.
%>
%> For every edge of a triangle the entries of g.markSideE0T indicate
%> whether the element is the first or second adjacent element with respect
%> to the edge view. @f$[K \times 3 \times 2]@f$
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
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
%> @author Balthasar Reuter, 2017
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