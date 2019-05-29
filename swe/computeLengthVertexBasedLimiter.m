% Computes a vector with correction factors for all elements to restrict
% a discrete vector field with two components, within the bounds of its
% euclidean norm in the surrounding centroids.
%
%===============================================================================
%> @file
%>
%> @brief Computes a vector with correction factors for all elements to
%>        restrict discrete vector field with two components, within the
%>        bounds of its euclidean norm in the surrounding centroids.
%===============================================================================
%>
%> @brief Computes a vector with correction factors for all elements to
%>        restrict discrete vector field with two components, within the
%>        bounds of its euclidean norm in the surrounding centroids.
%> 
%> This function is used as part of the vertex based slope limiter (cf. 
%> <code>applyLengthVertexBasedLimiter()</code>) to determine the
%> correction factors @f$\alpha_{ke}@f$ for all elements @f$T_k@f$.
%> 
%> It computes the bounds @f$\mathbf{u}_{ki}^\mathrm{max}@f$ 
%> (cf. <code>computeMinMaxV0TElementPatch()</code>)
%> using the given centroid values 
%> @f$\mathbf{u}_{kc} = \mathbf{u}_h(\mathbf{x}_{kc})@f$
%> and uses these bounds to obtain the correction factors by solving the
%> quadratic equation arising from the case that the inequality
%> @f[
%> || \mathbf{u}_{kc} + \alpha_{ke} (\mathbf{u}_{ki} - \mathbf{u}_{kc}) || 
%> \leq \mathbf{u}_{ki}^\mathrm{max}
%> @f]
%> is violated.
%>
%> @param  g            The lists describing the geometric and topological 
%>                      properties of a triangulation (see 
%>                      <code>generateGridData()</code>) 
%>                      @f$[1 \times 1 \text{ struct}]@f$
%> @param  valCentroid  The centroid values @f$\mathbf{u}_{kc}@f$. 
%>                      @f$[2 \times 1 \mathrm{cell}]@f$
%> @param  valV0T       The field values @f$\mathbf{u}_{ki}@f$ in the vertices. 
%>                      @f$[2 \times 1 \mathrm{cell}]@f$
%> @param  markV0TbdrD  <code>logical</code> arrays that mark each triangles
%>                      (Dirichlet boundary) vertices on which additional
%>                      function values should be considered during the 
%>                      slope limiting routine. @f$[K \times 3]@f$
%> @param  dataV0T      The vector field values for (Dirichlet boundary) vertices
%>                      specified by <code>markV0TbdrD</code>.
%>                      @f$[2 \times 1 \mathrm{cell}]@f$
%> @retval  alphaE      The correction factors for all elements. 
%>                      @f$[K\times 3]@f$
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
function alphaE = computeLengthVertexBasedLimiter(g, valCentroid, valV0T, markV0TbdrD, dataV0T)
% Check function arguments that are directly used
validateattributes(valCentroid, {'cell'}, {'size', [2 1]}, mfilename, 'valCentroid');
validateattributes(valV0T, {'cell'}, {'size', [2 1]}, mfilename, 'valV0T');
validateattributes(dataV0T, {'cell'}, {'size', [2 1]}, mfilename, 'dataV0T');

normCentroid = valCentroid{1}.^2+valCentroid{2}.^2;
minMaxV0T = computeMinMaxV0TElementPatch(g, normCentroid, markV0TbdrD, dataV0T{1}.^2 + dataV0T{2}.^2);

alphaEV0T = ones(g.numT,3);

tol = 1e-12;

for i = 1 : 3
  ind = valV0T{1}(:,i).^2 + valV0T{2}(:,i).^2 - minMaxV0T{2}(:,i) > tol;
  diff = { valV0T{1}(ind,i)-valCentroid{1}(ind), valV0T{2}(ind,i)-valCentroid{2}(ind) };
  normDiff = (diff{1}.^2 + diff{2}.^2);
  prod = - ( valCentroid{1}(ind) .* diff{1} + valCentroid{2}(ind) .* diff{2} );
  root = sqrt( prod.^2 - normDiff .* (normCentroid(ind) - minMaxV0T{2}(ind,i)) + tol );
  
  % Compute limiter parameter for each vertex in each triangle and ensure
  % that the value stays in [0,1]
  alphaEV0T(ind,i) = max(0, min(1, (prod + root) ./ normDiff));
end % for

% Each elements limiter parameter is the minimum of the three vertex values
alphaE = min(alphaEV0T,[],2);
end % function
