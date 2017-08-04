% Computes various additional fields needed throughout the advection
% problem presented in @ref RAWFK2016 . These fields only contain
% information based on the underlying grid and are therefore constant in
% time but used in every step. They are saved as part of the grid.

%===============================================================================
%> @file advection_implicit/computeDerivedGridData.m
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
%> g.areaE0TbdrD provides the elementwise products of the vectors 
%> <code>g.areaE0T</code>, and <code>g.markE0TbdrD</code>.
%> @f$[3 \text{ cell}]@f$
%> 
%> g.areaE0TbdrN provides the elementwise products of the vectors 
%> <code>g.areaE0T</code>, and <code>g.markE0TbdrN</code>.
%> @f$[K \times 3]@f$
%>
%> g.areaE0TbdrNotN provides the elementwise products of the vectors 
%> <code>g.areaE0T</code>, and <code>~g.markE0TbdrN</code>.
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
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
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
function g = computeDerivedGridData(g)
g.areaE0TbdrD = cell(3, 1);
g.areaE0TbdrN = g.markE0TbdrN .* g.areaE0T;
g.areaE0TbdrNotN = cell(3, 1);
for n = 1 : 3
  g.areaE0TbdrD{n} = g.areaE0T(:,n) .* g.markE0TbdrD(:, n);
  g.areaE0TbdrNotN{n} = g.areaE0T(:,n) .* ~g.markE0TbdrN(:, n);
end % for
end % function