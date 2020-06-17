% Computes a vector with correction factors for all elements to restrict
% a discrete function, given as vertex values, within the bounds of 
% surrounding centroid values.
%
%===============================================================================
%> @file
%>
%> @brief Computes a vector with correction factors for all elements to restrict
%>        a discrete function, given as vertex values, within the bounds of 
%>        surrounding centroid values.
%===============================================================================
%>
%> @brief Computes a vector with correction factors @f$\alpha_{ke}@f$ for 
%>        all elements to restrict a discrete function @f$c_h(\mathbf{x})@f$,
%>        given as vertex values, within the bounds of surrounding centroid values.
%> 
%> This function is used as part of vertex based slope limiters (cf. 
%> <code>applySlopeLimiterTaylorLinear()</code>) to determine the
%> correction factors @f$\alpha_{ke}@f$ for all elements @f$T_k@f$.
%> 
%> It computes the bounds @f$c_{ki}^\mathrm{min}@f$, 
%> @f$c_{ki}^\mathrm{max}@f$ (cf. <code>computeMinMaxV0TElementPatch()</code>)
%> using the given centroid values @f$c_{kc} = c_h(\mathbf{x}_{kc})@f$
%> and uses these bounds to obtain the correction factors by
%> @f[
%> \forall T_k\in\mathcal{T}_h\,,\quad  \alpha_{ke} :=
%>  \min_{i\,\in\, \{1, 2, 3\}} \left\{
%>  \begin{array}{llr}
%>    (c_{ki}^\mathrm{max} - c_{k\mathrm{c}})\big/(c_{ki} - c_{k\mathrm{c}})&\text{if}\;&c_{ki} > c_{ki}^\mathrm{max}\\
%>     1                  &\text{if}\;&c_{ki}^\mathrm{min}\le c_{ki}\le c_{ki}^\mathrm{max}\\
%>    (c_{ki}^\mathrm{min} - c_{k\mathrm{c}})\big/(c_{ki} - c_{k\mathrm{c}})&\text{if}\;&c_{ki} < c_{ki}^\mathrm{min}
%>  \end{array}\right\}\;,
%> @f]
%> where @f$c_{ki} = c_h(\mathbf{x}_{kc})@f$ are the function values on the
%> vertices of triangle @f$T_k@f$.
%>
%> @param  g           The lists describing the geometric and topological 
%>                     properties of a triangulation (see 
%>                     <code>generateGridData()</code>) 
%>                     @f$[1 \times 1 \text{ struct}]@f$
%> @param  valCentroid The centroid values @f$c_{kc}@f$. @f$[K\times 1]@f$
%> @param  valV0T      The function values @f$c_{ki}@f$. @f$[K\times 3]@f$
%> @param  markV0TbdrD <code>logical</code> array that marks each triangles
%>                     (Dirichlet boundary) vertices on which the values given
%>                     in dataV0T should be included in the computation of
%>                     minimum/maximum @f$[K \times 3]@f$
%> @param  dataV0T     A list of values to be included in the determination
%>                     of minimum/maximum. @f$[K\times 3]@f$
%> @retval  alphaE     The correction factors for all elements. 
%>                     @f$[K\times 3]@f$
%> @retval  minMaxV0T  Two matrices with minimum or maximum centroid values,
%>                     respectively, of the patch of elements surrounding each
%>                     vertex of each element as computed by 
%>                     <code>computeMinMaxV0TElementPatch()</code>
%>                     @f$[2 \times 1 \mathrm{cell}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Florian Frank, Balthasar Reuter, Vadym Aizinger
%>                      Modified 09/02/16 by Hennes Hajduk
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
function [alphaE, minMaxV0T] = computeVertexBasedLimiter(g, valCentroid, valV0T, markV0TbdrD, dataV0T)
% Check function arguments that are directly used
n_edges = size(g.V0T, 2);
validateattributes(valCentroid, {'numeric'}, {'size', [g.numT 1]}, mfilename, 'valCentroid');
validateattributes(dataV0T, {'numeric'}, {'size', [g.numT n_edges]}, mfilename, 'dataV0T');

% Determine for each vertex the min- and max-values of the centroids in 
% elements adjacent to the vertex
minMaxV0T = computeMinMaxV0TElementPatch(g, valCentroid, markV0TbdrD, dataV0T);

% Compute deviance from centroid values for vertex- and min/max-values
diffV0TCentroid = valV0T - repmat(valCentroid, [1 n_edges]);
diffMinCentroid = minMaxV0T{1} - repmat(valCentroid, [1 n_edges]);
diffMaxCentroid = minMaxV0T{2} - repmat(valCentroid, [1 n_edges]);

% Find out if positive or negative deviance from centroid value
tol = 1.e-8;
markNeg = diffV0TCentroid < diffMinCentroid + tol;
markPos = diffV0TCentroid > diffMaxCentroid - tol;

% Compute limiter parameter for each vertex in each triangle and ensure
% that the value stays in [0,1]
alphaEV0T = ones(g.numT, n_edges);
alphaEV0T(markNeg) = max(0, min(1, diffMinCentroid(markNeg) ./ (diffV0TCentroid(markNeg) - tol) ) );
alphaEV0T(markPos) = max(0, min(1, diffMaxCentroid(markPos) ./ (diffV0TCentroid(markPos) + tol) ) );

% Each elements limiter parameter is the minimum of the three vertex values
alphaE = min(alphaEV0T,[],2);
end % function