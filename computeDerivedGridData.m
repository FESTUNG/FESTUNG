% Computes various additional fields needed throughout the diffusion
% problem. These fields only contain information based on the underlying
% grid and are therefore constant in time but used in every step. They are 
% saved as part of the grid.

%===============================================================================
%> @file diffusion/computeDerivedGridData.m
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
%> g.areaNuE0TE0T provides the elementwise products of the matrices 
%> <code>g.markE0TE0T</code> with the vectors of the elementwise product of
%> <code>g.areaE0T</code>, and <code>g.nuE0T</code>.
%> @f$[3 \times 3 \times 2 \text{ cell}]@f$
%>
%> g.areaNuE0T provides the elementwise products of the vectors 
%> <code>g.areaE0T</code>, and <code>g.nuE0T</code>.
%> @f$[3 \times 2 \text{ cell}]@f$
%>
%> g.areaNuE0Tint provides the elementwise products of the vectors 
%> <code>g.areaE0T</code>, <code>g.nuE0T</code>, and <code>markE0Tint</code>.
%> @f$[3 \times 2 \text{ cell}]@f$
%>
%> g.areaNuE0TbdrD provides the elementwise products of the vectors 
%> <code>g.areaE0T</code>, <code>g.nuE0T</code>, and <code>markE0TbdrD</code>.
%> @f$[3 \times 2 \text{ cell}]@f$
%>
%> g.areaNuE0TbdrD provides the elementwise products of the vectors 
%> <code>g.areaE0T</code>, <code>g.nuE0T</code>, and <code>markE0TbdrN</code>.
%> @f$[3 \times 2 \text{ cell}]@f$
%>
%> g.areaE0TbdrN provides the elementwise products of the vectors 
%> <code>g.areaE0T</code>, and <code>markE0TbdrN</code>.
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
g.areaNuE0TE0T = cell(3,3,2);
g.areaNuE0T = cell(3,2);
g.areaNuE0Tint = cell(3,2);
g.areaNuE0TbdrD = cell(3,2);
g.areaNuE0TbdrN = cell(3,2);
g.areaE0TbdrN = cell(3,1);
for nn = 1 : 3
  for np = 1 : 3
    for m = 1 : 2
      g.areaNuE0TE0T{nn,np,m} = bsxfun(@times, g.markE0TE0T{nn,np}, g.areaE0T(:,nn).*g.nuE0T(:,nn,m));
    end % for
  end % for
  for m = 1 : 2
    g.areaNuE0T{nn,m} = g.areaE0T(:,nn).*g.nuE0T(:,nn,m);
    g.areaNuE0Tint{nn,m} = g.areaNuE0T{nn,m} .* g.markE0Tint(:, nn);
    g.areaNuE0TbdrD{nn,m} = g.areaNuE0T{nn,m} .* g.markE0TbdrD(:,nn);
    g.areaNuE0TbdrN{nn,m} = g.areaNuE0T{nn,m} .* g.markE0TbdrN(:,nn);
  end % for
  g.areaE0TbdrN{nn} = g.markE0TbdrN(:, nn) .* g.areaE0T(:,nn);
end % for
end % function