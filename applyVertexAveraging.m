% Takes a DG function and smoothens it by averaging the vertex values from 
% different triangles.

%===============================================================================
%> @file applyVertexAveraging.m
%>
%> @brief Takes a DG function and smoothens it by averaging the vertex values 
%>        from different triangles.
%===============================================================================
%>
%> @brief Takes a DG function and smoothens it by averaging the vertex values 
%>        from different triangles.
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%>
%> @param  data       A representation of the discrete function before the smoothening
%>                    @f$d_h(\mathbf(x))@f$, e.g., as computed by 
%>                    <code>projectFuncCont2DataDisc()</code>
%>                    @f$[K \times N]@f$
%>
%> @param  averagingOperator matrix that is non-zero in row i, column 
%>					  3(k-1)+j (1 <= k <= K, 1 <= j <= 3) iff global vertex i is the j-th 
%>                    local vertex of element k. The sum of each row is one. 
%>                    @f$[numV \times 3K]@f$
%>
%> @param phiAtVertex The values of the first three basis functions in each of 
%>                    the vertices on the reference element.
%>					  @f$[3 \times 3]@f$
%>
%> @retval data       A representation of the discrete function after the smoothening
%>                    @f$d_h(\mathbf(x))@f$, e.g., as computed by 
%>                    <code>projectFuncCont2DataDisc()</code>
%>                    @f$[K \times N]@f$
%> @retval dataV      A representation of the discrete function after the smoothening
%>                    @f$d_h(\mathbf(x))@f$ given by the values in the grid vertices
%>                    @f$[numV \times 1]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Florian Frank, Balthasar Reuter, Vadym Aizinger
%>
%> @author Hennes Hajduk, 2017
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
function [data, dataV] = applyVertexAveraging(g, data, averagingOperator, phiAtVertex)
dataLagr = projectDataDisc2DataLagr(data, 1);
dataV = averagingOperator * reshape(dataLagr.', 3*g.numT, 1);
dataV0T = dataV(g.V0T);
data = dataV0T / phiAtVertex;
end % function