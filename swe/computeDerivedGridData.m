% Computes various additional fields needed repeatedly throughout the 
% shallow water problem as it is solved in @ref @ref HHAR2018 . These 
% fields only contain information based on the specific grid and are 
% therefore constant in time. They are stored as part of the grid.

%===============================================================================
%> @file
%>
%> @brief Computes various additional fields needed repeatedly throughout
%>				the shallow water problem as it is solved in @ref @ref HHAR2018.
%> 				These fields only contain information based on the specific grid
%> 				and are therefore constant in time. They are stored as part of 
%>        the grid.
%===============================================================================
%>
%> @brief Computes various additional fields needed repeatedly throughout
%>				the shallow water problem as it is solved in @ref @ref HHAR2018.
%> 				These fields only contain information based on the specific grid
%> 				and are therefore constant in time. They are stored as part of 
%>        the grid.
%>
%> g.markE0T provides the vector equal to <code>g.markE0TE0T</code> 
%> summed over each row.
%> @f$[3 \times 3 \text{ cell}]@f$
%> 
%> For every vertex of a triangle the columns of g.markV0TT0V indicate which triangle 
%> corresponds to this vertex.
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @retval g          The lists describing the geometric and topological 
%>                    properties of a triangulation enriched with the 
%>                    precomputed fields.
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2018 Florian Frank, Balthasar Reuter, Vadym Aizinger
%> 
%> @author Hennes Hajduk, 2018
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
g.markE0T = cell(3,3);
g.markV0TT0V = cell(3,1	);
for nn = 1 : 3
	for np = 1 : 3
    g.markE0T{nn,np} = g.markE0TE0T{nn,np} * ones(g.numT,1);
	end % for
  g.markV0TT0V{nn} = g.markV0TV0T{nn, 1} | g.markV0TV0T{nn, 2} | g.markV0TV0T{nn, 3};
end % for
end % function
