% Computes various additional fields needed throughout the diffusion
% problem. These fields only contain information based on the underlying
% grid and are therefore constant in time but used in every step. They are 
% saved as part of the grid.

%===============================================================================
%> @file
%>
%> @brief Computes various additional fields needed throughout the diffusion
%>        problem presented in @ref FRAK2015 . These fields only contain
%>        information based on the underlying grid and are therefore constant in
%>        time but used in every step. They are saved as part of the grid.
%===============================================================================
%>
%> @brief Computes various additional fields needed throughout the diffusion
%>        problem presented in @ref FRAK2015 . These fields only contain
%>        information based on the underlying grid and are therefore constant in
%>        time but used in every step. They are saved as part of the grid.
%>
%> g.areaNuE0T provides the elementwise products of the vectors 
%> <code>g.areaE0T</code>, and <code>g.nuE0T</code>.
%> @f$[K \times 3 \times 2]@f$
%>
%> g.areaNuE0Tint provides the elementwise products of the vectors 
%> <code>g.areaE0T</code>, <code>g.nuE0T</code>, and <code>markE0Tint</code>.
%> @f$[2 \times 1 \text{ cell}]@f$
%>
%> g.areaNuE0TbdrD provides the elementwise products of the vectors 
%> <code>g.areaE0T</code>, <code>g.nuE0T</code>, and <code>markE0TbdrD</code>.
%> @f$[2 \times 1 \text{ cell}]@f$
%>
%> g.areaNuE0TbdrD provides the elementwise products of the vectors 
%> <code>g.areaE0T</code>, <code>g.nuE0T</code>, and <code>markE0TbdrN</code>.
%> @f$[2 \times 1 \text{ cell}]@f$
%>
%> g.areaE0TbdrN provides the elementwise products of the vectors 
%> <code>g.areaE0T</code>, and <code>markE0TbdrN</code>.
%> @f$[K \times 3]@f$
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
%> @copyright 2014-2015 Florian Frank, Balthasar Reuter, Vadym Aizinger
%>
%> @author Hennes Hajduk, 2016
%> @author Balthasar Reuter, 2018
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
K = g.numT;
g.areaNuE0T = zeros(K,3,2);
g.areaNuE0Tint = cell(2,1);
g.areaNuE0TbdrD = cell(2,1);
g.areaNuE0TbdrN = cell(2,1);
g.areaE0TbdrN = g.markE0TbdrN .* g.areaE0T;
for m = 1 : 2
  g.areaNuE0T(:, :, m) = g.areaE0T .* g.nuE0T(:, :, m);
  g.areaNuE0Tint{m} = g.markE0Tint .* g.areaNuE0T(:, :, m);
  g.areaNuE0TbdrD{m} = g.markE0TbdrD .* g.areaNuE0T(:, :, m);
  g.areaNuE0TbdrN{m} = g.markE0TbdrN .* g.areaNuE0T(:, :, m);
end % for m
end % function