% Computes various additional fields needed throughout the Shallow-Water
% problem presented in ... . These fields only contain
% information based on the underlying grid and are therefore constant in
% time but used in every step. They are stored as part of the grid.

%===============================================================================
%> @file computeDerivedGridData.m
%>
%> @brief Computes various additional fields needed throughout the Shallow-Water
%>				problem presented in ... . These fields only contain
%>				information based on the underlying grid and are therefore constant in
%>				time but used in every step. They are stored as part of the grid.
%===============================================================================
%>
%> @brief Computes various additional fields needed throughout the Shallow-Water
%>				problem presented in ... . These fields only contain
%>				information based on the underlying grid and are therefore constant in
%>				time but used in every step. They are stored as part of the grid.
%>
%> g.areaE0Tint provides the elementwise products of the vectors 
%> <code>g.areaE0T</code>, and <code>markE0Tint</code>.
%> @f$[3 \times 1 \text{ cell}]@f$
%>
%> g.areaE0TbdrL provides the elementwise products of the vectors 
%> <code>g.areaE0T</code>, and <code>markE0TbdrL</code>.
%> @f$[3 \times 1 \text{ cell}]@f$
%>
%> g.areaE0TbdrRI provides the elementwise products of the vectors 
%> <code>g.areaE0T</code>, and <code>markE0TbdrRI</code>.
%> @f$[3 \times 1 \text{ cell}]@f$
%>
%> g.areaE0TbdrOS provides the elementwise products of the vectors 
%> <code>g.areaE0T</code>, and <code>markE0TbdrOS</code>.
%> @f$[3 \times 1 \text{ cell}]@f$
%>
%> g.areaNuE0Tint provides the elementwise products of the vectors 
%> <code>g.areaE0T</code>, <code>g.nuE0T</code> and <code>markE0Tint</code>.
%> @f$[3 \times 2 \text{ cell}]@f$
%>
%> g.areaNuE0TbdrL provides the elementwise products of the vectors 
%> <code>g.areaE0T</code>, <code>g.nuE0T</code> and <code>markE0TbdrL</code>.
%> @f$[3 \times 2 \text{ cell}]@f$
%>
%> g.areaNuE0TbdrRA provides the elementwise products of the vectors 
%> <code>g.areaE0T</code>, <code>g.nuE0T</code> and <code>markE0TbdrRA</code>.
%> @f$[3 \times 2 \text{ cell}]@f$
%>
%> g.areaNuE0TbdrRI provides the elementwise products of the vectors 
%> <code>g.areaE0T</code>, <code>g.nuE0T</code> and <code>markE0TbdrRI</code>.
%> @f$[3 \times 2 \text{ cell}]@f$
%>
%> g.areaNuE0TbdrOS provides the elementwise products of the vectors 
%> <code>g.areaE0T</code>, <code>g.nuE0T</code> and <code>markE0TbdrOS</code>.
%> @f$[3 \times 2 \text{ cell}]@f$
%>
%> g.markE0T provides the vector equal to <code>g.markE0TE0T</code> 
%> summed over each row.
%> @f$[3 \times 3 \text{ cell}]@f$
%>
%> g.areaMarkE0T provides the elementwise products of the vectors 
%> <code>g.areaE0T</code>, and the vector equal to <code>g.markE0TE0T</code> 
%> summed over each row.
%> @f$[3 \times 3 \text{ cell}]@f$
%>
%> g.areaMarkE0TE0T provides the elementwise products of each column of the
%> matrices <code>g.markE0TE0T</code> with the vectors that make up the 
%> three columns of <code>g.areaE0T</code>. @f$[3 \times 3 \text{ cell}]@f$
%>
%> g.areaNuMarkE0TE0T provides the elementwise products of each column of the
%> matrices <code>g.markE0TE0T</code> with the vectors of the elementwise product
%> of the columns of <code>g.areaE0T</code>, and <code>g.nuE0T</code> for each 
%> normal component. @f$[3 \times 3 \times 2 \text{ cell}]@f$
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
%> @copyright 2014-2016 Hennes Hajduk, Florian Frank, Balthasar Reuter, Vadym Aizinger
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
g.areaE0Tint				= cell(3,1	);
g.areaE0TbdrL 			= cell(3,1	);
g.areaE0TbdrRI			= cell(3,1	);
g.areaE0TbdrOS			= cell(3,1	);
g.areaNuE0Tint			= cell(3,2	);
g.areaNuE0TbdrL			= cell(3,2	);
g.areaNuE0TbdrRA		= cell(3,2  );
g.areaNuE0TbdrRI		= cell(3,2  );
g.areaNuE0TbdrOS		= cell(3,2	);
g.areaMarkE0T				= cell(3,3	);
g.areaMarkE0TE0T		= cell(3,3	);
g.areaNuMarkE0TE0T	= cell(3,3,2);
areaNuE0T						= cell(3,2  );
for nn = 1 : 3
	g.areaE0Tint	{nn} = g.areaE0T(:,nn) .* g.markE0Tint	(:,nn);
	g.areaE0TbdrL {nn} = g.areaE0T(:,nn) .* g.markE0TbdrL (:,nn);
  g.areaE0TbdrRI{nn} = g.areaE0T(:,nn) .* g.markE0TbdrRI(:,nn);
	g.areaE0TbdrOS{nn} = g.areaE0T(:,nn) .* g.markE0TbdrOS(:,nn);
	for m = 1 : 2
    areaNuE0T{nn,m} = g.areaE0T(:,nn) .* g.nuE0T(:,nn,m);
		g.areaNuE0Tint	{nn,m} = areaNuE0T{nn,m} .* g.markE0Tint	(:,nn);
		g.areaNuE0TbdrL {nn,m} = areaNuE0T{nn,m} .* g.markE0TbdrL	(:,nn);
		g.areaNuE0TbdrRA{nn,m} = areaNuE0T{nn,m} .* g.markE0TbdrRA(:,nn);
		g.areaNuE0TbdrRI{nn,m} = areaNuE0T{nn,m} .* g.markE0TbdrRI(:,nn);
		g.areaNuE0TbdrOS{nn,m} = areaNuE0T{nn,m} .* g.markE0TbdrOS(:,nn);
	end % for
	for np = 1 : 3
    g.markE0T		    {nn,np} = g.markE0TE0T{nn,np} * ones(g.numT,1);
		g.areaMarkE0T		{nn,np} = g.areaE0T(:,nn) .* g.markE0T{nn,np};
		g.areaMarkE0TE0T{nn,np} = bsxfun( @times, g.markE0TE0T{nn, np}, g.areaE0T(:,nn) );
		for m = 1 : 2
			g.areaNuMarkE0TE0T{nn,np,m} = bsxfun(@times, g.markE0TE0T{nn, np}, areaNuE0T{nn,m});
		end % for
	end % for
end % for
end % function
