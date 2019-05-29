% Computes a vector with correction factors for all elements to restrict
% a two element vector of discrete functions within the bounds of 
% surrounding centroid values and a corresponding rotation matrix.
%
%===============================================================================
%> @file
%>
%> @brief Computes a vector with correction factors for all elements to restrict
%>        a two element vector of discrete functions within the bounds of 
%>        surrounding centroid values and a corresponding rotation matrix.
%===============================================================================
%>
%> @brief Computes a vector with correction factors for all elements to restrict
%>        a two element vector of discrete functions within the bounds of 
%>        surrounding centroid values and a corresponding rotation matrix.
%> 
%> This function is used as part of vertex based slope limiters (cf. 
%> <code>applySlopeLimiterTaylorLinear()</code>) to determine the
%> correction factors @f$\alpha_{ke}@f$ for all elements @f$T_k@f$.
%> 
%> It first computes a rotation matrix, according to a specified method and
%> then the minima and maxima of 
%> @f$\mathbf{q}_j^k \cdot (\mathbf{u}_{l\mathrm{c}} - \mathbf{u}_{k\mathrm{c}})@f$
%> where @f$l@f$ is the index of an element sharing a node with element 
%> @f$T_k$ for each node of each element
%> (cf. <code>computeMinMaxVectorV0TElementPatch()</code>).
%> The limiting directions @f$\mathbf{q}_j^k@f$ given by the rows of the 
%> rotaion matrix, and the centroid values are @f$c_{kc} = c_h(\mathbf{x}_{kc})@f$.
%> Then, the correction factors are obtained by
%> @f[
%> \forall T_k\in\mathcal{T}_h\,,\quad  \alpha_{ke} :=
%>  \min_{i\,\in\, \{1, 2, 3\}} \left\{
%>  \begin{array}{llr}
%>    (\mathrm{max} \mathbf{q}_j^k \cdot (\mathbf{u}_{l\mathrm{c}} - \mathbf{u}_{k\mathrm{c}}))\big/
%>      (\mathbf{q}_j^k \cdot (\mathbf{u}_{ki} - \mathbf{u}_{k\mathrm{c}}))
%>      &\text{if}\;&(\mathbf{q}_j^k \cdot (\mathbf{u}_{ki} - \mathbf{u}_{k\mathrm{c}})) > 0\\
%>     1       &\text{if}\;&(\mathbf{q}_j^k \cdot (\mathbf{u}_{ki} - \mathbf{u}_{k\mathrm{c}})) = 0\\
%>    (\mathrm{min} \mathbf{q}_j^k \cdot (\mathbf{u}_{l\mathrm{c}} - \mathbf{u}_{k\mathrm{c}}))\big/
%>      (\mathbf{q}_j^k \cdot (\mathbf{u}_{ki} - \mathbf{u}_{k\mathrm{c}}))
%>      &\text{if}\;&(\mathbf{q}_j^k \cdot (\mathbf{u}_{ki} - \mathbf{u}_{k\mathrm{c}})) < 0
%>  \end{array}\right\}\;,
%> @f]
%> where @f$\mathbf{u}_{ki} = \mathbf{u}_h(\mathbf{x}_{ki})@f$ are the function values on the
%> vertices of triangle @f$T_k@f$.
%> The rotation matrix is computed in <code>computeRotationMatrix()</code>.
%>
%> @param  g            The lists describing the geometric and topological 
%>                      properties of a triangulation (see 
%>                      <code>generateGridData()</code>) 
%>                      @f$[1 \times 1 \text{ struct}]@f$
%> @param  dataTaylor   The representation matrices of the unlimited vector field 
%>                      @f$\mathbf{u}_h@f$. @f$[2 \times 1 \mathrm{cell}]@f$
%> @param  valV0T       The function values @f$\mathbf{u}_{ki}@f$. 
%>                      @f$[2 \times 1 \mathrm{cell}]@f$
%> @param  markV0TbdrD  <code>logical</code> array that marks each triangles
%>                      (Dirichlet boundary) vertices on which the values given
%>                      in dataV0T should be included in the computation of
%>                      minimum/maximum @f$[K \times 3]@f$
%> @param  dataV0T      A list of values to be included in the determination
%>                      of minimum/maximum. @f$[2 \times 1 \mathrm{cell}]@f$
%> @param  variant      The variant that is used to compute rotation
%>                      matrices for directional vector limiting.
%> @param  angle        The rotation angle to be used in case a fixed
%>                      transformation is desired to compute the limiting 
%>                      directions.
%> @retval  alphaE      The correction factors for all elements for each
%>                      limiting direction
%>                      @f$ @f$[2 \times 1 \mathrm{cell}]@f$.
%> @retval Q            The rotation matrixces for each element.
%>                      @f$[2 \times 2 \mathrm{cell}]@f$
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
function [alphaE, Q] = computeVectorVertexBasedLimiter(g, dataTaylor, valV0T, markV0TbdrD, dataV0T, variant, angle)

Q = computeRotationMatrix(g, {dataTaylor{1}(:,2), dataTaylor{1}(:,3); dataTaylor{2}(:,2), dataTaylor{2}(:,3)}, variant, angle);

valCentroid = {dataTaylor{1}(:,1); dataTaylor{2}(:,1)};

alphaE = cell(2,1);

for i = 1:2
	minMaxV0T = computeMinMaxVectorV0TElementPatch(g, Q(i,:), valCentroid, markV0TbdrD, dataV0T);
	
	diff = bsxfun(@times, valV0T{1} - repmat(valCentroid{1}, [1 3]), Q{i,1}) + ...
				 bsxfun(@times, valV0T{2} - repmat(valCentroid{2}, [1 3]), Q{i,2});
	
	% Find out if positive or negative deviance from centroid value
	tol = 1.e-8;
	markNeg = diff < minMaxV0T{1} + tol;
	markPos = diff > minMaxV0T{2} - tol;
	
	% Compute limiter parameter for each vertex in each triangle and ensure
	% that the value stays in [0,1]
	alphaEV0T = ones(g.numT,3);
	alphaEV0T(markNeg) = max(0, min(1, minMaxV0T{1}(markNeg) ./ (diff(markNeg) - tol) ) );
	alphaEV0T(markPos) = max(0, min(1, minMaxV0T{2}(markPos) ./ (diff(markPos) + tol) ) );
	
	% Each elements limiter parameter is the minimum of the three vertex values
	alphaE{i} = min(alphaEV0T,[],2);
end % for

end % function
