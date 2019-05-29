% Applies a slope limiter for the euclidean norm to a two component vector of
% functions.
%===============================================================================
%> @file
%>
%> @brief Applies a slope limiter for the euclidean norm to a two component
%>        vector of functions.
%===============================================================================
%>
%> @brief Applies a slope limiter for the euclidean norm to a two component
%>        vector of functions.
%>
%> It applies the vertex based limiter to a vector field in order to satisfy
%> maximum principles for its length.
%>
%> @param  g            The lists describing the geometric and topological 
%>                      properties of a triangulation (see 
%>                      <code>generateGridData()</code>) 
%>                      @f$[1 \times 1 \text{ struct}]@f$
%> @param  dataTaylor   The representation matrix of the unlimited vector 
%>                      field @f$\mathbf{u}_h@f$. 
%>                      @f$[1 \times 2 \mathrm{cell}]@f$
%> @param  valV0T       The function values @f$\mathbf{u}_{ki}@f$. 
%>                      @f$[2 \times 1 \mathrm{cell}]@f$
%> @param  markV0TbdrD  <code>logical</code> arrays that mark each triangles
%>                      (Dirichlet boundary) vertices on which additional
%>                      function values should be considered during the 
%>                      slope limiting routine. @f$[K \times 3]@f$
%> @param  dataV0T      The vector field values for (Dirichlet boundary) vertices
%>                      specified by <code>markV0TbdrD</code>. @f$[K \times 3]@f$
%> @retval dataTaylor   The representation matrix of the limited vector 
%>                      field @f$\mathbf{u}_h@f$. 
%>                      @f$[1 \times 2 \mathrm{cell}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2018 Florian Frank, Balthasar Reuter, Vadym Aizinger
%>
%> @author Hennes Hajduk, 2018.
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
function dataTaylor = applyLengthVertexBasedLimiter(g, dataTaylor, valV0T, markV0TbdrD, dataV0T)
% Check function arguments that are directly used
validateattributes(dataTaylor, {'cell'}, {'size', [1 2]}, mfilename, 'dataTaylor');

%% Limit first order derivative terms
% Compute limiter parameter for each vertex
alphaE = computeLengthVertexBasedLimiter(g, {dataTaylor{1}(:,1); dataTaylor{2}(:,1)}, valV0T, markV0TbdrD, dataV0T);

% Apply limiter to first order terms, set all higher order terms to zero in the case of limiting
dataTaylor{1} = [dataTaylor{1}(:,1), bsxfun(@times, alphaE, dataTaylor{1}(:, 2:3)), bsxfun(@times, alphaE == 1, dataTaylor{1}(:, 4:end))];
dataTaylor{2} = [dataTaylor{2}(:,1), bsxfun(@times, alphaE, dataTaylor{2}(:, 2:3)), bsxfun(@times, alphaE == 1, dataTaylor{2}(:, 4:end))];
end % function
