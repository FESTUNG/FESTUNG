% Computes various additional fields needed throughout the Shallow-Water
% problem presented in ... . These fields only contain
% information based on the underlying grid and are therefore constant in
% time but used in every step. They are saved as part of the grid.
%
%===============================================================================
%> @file computeDerivedGridDataSWE.m
%>
%> @brief Computes various additional fields needed throughout the Shallow-Water
%>				problem presented in ... . These fields only contain
%>				information based on the underlying grid and are therefore constant in
%>				time but used in every step. They are saved as part of the grid.
%===============================================================================
%>
%> @brief Computes various additional fields needed throughout the Shallow-Water
%>				problem presented in ... . These fields only contain
%>				information based on the underlying grid and are therefore constant in
%>				time but used in every step. They are saved as part of the grid.
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
%> g.areaNuE0TbdrOS provides the elementwise products of the vectors 
%> <code>g.areaE0T</code>, <code>g.nuE0T</code> and <code>markE0TbdrOS</code>.
%> @f$[3 \times 2 \text{ cell}]@f$
%>
%> g.areaE0TsumCols provides the elementwise products of the vectors 
%> <code>g.areaE0T</code>, and the vector equal to <code>g.markE0TE0T</code> 
%> summed over each row.
%> @f$[3 \times 3 \text{ cell}]@f$
%>
%> g.areaNuE0TE0T provides the elementwise products of each column of the
%> matrices <code>g.markE0TE0T</code> with the vectors of the elementwise product
%> of <code>g.areaE0T</code>, and <code>g.nuE0T</code>.
%> @f$[3 \times 3 \times 2 \text{ cell}]@f$
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0Tint <code>logical</code> arrays that mark each triangles
%>                    (interior) edges on which the matrix blocks should be
%>                    assembled @f$[K \times 3]@f$
%> @param  markE0TbdrL <code>logical</code> arrays that mark each triangles
%>                    (land boundary) edges on which the vector entries should be
%>                    assembled @f$[K \times 3]@f$
%> @param  markE0TbdrOS <code>logical</code> arrays that mark each triangles
%>                    (open sea boundary) edges on which the vector entries should be
%>                    assembled @f$[K \times 3]@f$
%> @retval g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
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
function g = computeDerivedGridDataSWE(g, markE0Tint, markE0TbdrL, markE0TbdrRA, markE0TbdrRI, markE0TbdrOS)
g.areaE0TbdrL				= cell(3,1	);
g.areaE0TbdrOS			= cell(3,1	);
g.areaNuE0Tint			= cell(3,2	);
g.areaNuE0TbdrL			= cell(3,2	);
g.areaNuE0TbdrRA		= cell(3,2  );
g.areaNuE0TbdrRI		= cell(3,2  );
g.areaNuE0TbdrOS		= cell(3,2	);
g.areaE0TsumCols  	= cell(3,3	);
g.areaNuE0TE0T      = cell(3,3,2);
areaNuE0T						= cell(3,2  );
for nn = 1 : 3
	g.areaE0TbdrL {nn} = g.areaE0T(:,nn) .* markE0TbdrL (:,nn);
	g.areaE0TbdrOS{nn} = g.areaE0T(:,nn) .* markE0TbdrOS(:,nn);
	for m = 1 : 2
			areaNuE0T			{nn,m} = g.areaE0T(:,nn) .* g.nuE0T(:,nn,m);
		g.areaNuE0Tint	{nn,m} = areaNuE0T{nn,m} .* markE0Tint	(:,nn);
		g.areaNuE0TbdrL {nn,m} = areaNuE0T{nn,m} .* markE0TbdrL	(:,nn);
		g.areaNuE0TbdrRA{nn,m} = areaNuE0T{nn,m} .* markE0TbdrRA(:,nn);
		g.areaNuE0TbdrRI{nn,m} = areaNuE0T{nn,m} .* markE0TbdrRI(:,nn);
		g.areaNuE0TbdrOS{nn,m} = areaNuE0T{nn,m} .* markE0TbdrOS(:,nn);
	end % for
	for np = 1 : 3
		g.areaE0TsumCols{nn,np} = g.areaE0T(:,nn) .* (g.markE0TE0T{nn,np} * ones(g.numT,1));
		for m = 1 : 2
			g.areaNuE0TE0T{nn,np,m} = bsxfun(@times, g.markE0TE0T{nn, np}, areaNuE0T{nn,m});
		end % for
	end % for
end % for
end % function
