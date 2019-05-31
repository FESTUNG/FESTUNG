% Modifies all required grid data in order to form periodic boundary conditions.
%
%===============================================================================
%> @file
%>
%> @brief Modifies all required grid data in order to form periodic boundary 
%>        conditions based on a list of edges of the grid that shall be
%>        identified with each other and an optional list of identified nodes.
%===============================================================================
%>
%> @brief Modifies all required grid data in order to form periodic boundary 
%>        conditions based on a list of edges of the grid that shall be
%>        identified with each other and an optional list of identified nodes.
%>
%> It modifies the topological edge and vertex data, resulting in a grid
%> over a rectangular domain with periodic boundary conditions. 
%>
%> @param  g                    The lists describing the geometric and
%>                              topological properties of a triangulation
%>                              (see <code>generateGridData()</code>)
%>                              @f$[1 \times 1 \text{ struct}]@f$
%> @param  edges                A list of edges indices describing which
%>                              edges shall be identified. 
%>                              @f$[numPeriodicEdges \times N]@f$
%> @param  nodes                A vector of node indices that correspond to
%>                              the same grid point.
%> @retval isPeriodicMarkV0TV0T Logical that indicates if the field
%>                              markV0TV0T has already been computed.
%> @retval g                    The modified lists describing the geometric
%>                              and topological properties of a triangulation
%>                              (see <code>generateGridData()</code>) 
%>                              @f$[1 \times 1 \text{ struct}]@f$
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
function g = modifyGridDataPeriodic(g, edges, nodes, isPeriodicMarkV0TV0T)

[numPerE, l] = size(edges);
assert(l == 2, 'Edge identification requires pairs of edge indices.');

V0T = g.V0T;
for i = 1:length(nodes)
  V0T(V0T==nodes(i)) = 0;
end % for

ctr = -1;

for i = 1:numPerE
	e1 = edges(i,1);
	e2 = edges(i,2);
	if ~isequal(g.T0E([e1 e2],2), [0; 0]) || (g.idE(e1) == 0) || (g.idE(e2) == 0)
		error('Interior edges cannot be set to periodic boundary edges.')
	end
	if abs(g.areaE(e1) - g.areaE(e2)) > 1.e-10
		error('Corresponding edges have different lengths. Change position of vertices!')
	end
	if (g.nuE(e1,1) ~= -g.nuE(e2,1)) || (g.nuE(e1,1) ~= -g.nuE(e2,1))
		warning('Corresponding edges do not have corresponding normals.')
	end
	
	T0E = g.T0E([e1 e2], 1);
	g.T0E([e2 e1],2) = T0E;
	
	locE = g.E0E([e1 e2], 1);
	g.E0E([e2 e1],2) = locE;
	
	g.markE0TE0T{locE(1), locE(2)}(T0E(1), T0E(2)) = 1;
	g.markE0TE0T{locE(2), locE(1)}(T0E(2), T0E(1)) = 1;
	
	V0E = g.V0E([e1 e2],:);
	
	g.V2T(V0E(1,2), V0E(1,1)) = T0E(2);
	g.V2T(V0E(2,2), V0E(2,1)) = T0E(1);
  
	V0T(V0T==V0E(1,1)) = ctr;
	V0T(V0T==V0E(2,2)) = ctr;
	V0T(V0T==V0E(1,2)) = ctr-1;
	V0T(V0T==V0E(2,1)) = ctr-1;
	ctr = ctr-2;
  
  g.idE([e1 e2]) = 0;
end % for

g.idE0T = g.idE(g.E0T);

if ~isPeriodicMarkV0TV0T
  for nn = 1 : 3
    for np = 1 : 3
      g.markV0TV0T{nn,np} = sparse(bsxfun(@eq, V0T(:,nn), V0T(:,np)'));
    end % for
  end % for
end % if

end % function
