% Applies the hierarchical vertex based slope limiter to a discrete 
% function given in Taylor basis.

%===============================================================================
%> @file
%>
%> @brief Applies the hierarchical vertex based slope limiter to a discrete 
%>        function given in Taylor basis.
%===============================================================================
%>
%> @brief Applies the hierarchical vertex based slope limiter to a discrete 
%>        function given in Taylor basis.
%>
%> This is an extension of the linear limiter (cf.
%> <code>applySlopeLimiterTaylorLinear()</code>). 
%> Here, a correction factor @f$\alpha_{ke}^{(q)}@f$ is determined for and 
%> applied to all derivatives of order @f$q@f$.
%>
%> Let @f$\mathcal{A}_q:=\{\mathbf{a}\in\mathbb{N}_0^2\,\big|\,|\mathbf{a}|=q\}@f$ 
%> be the set of all two-dimensional multi-indices of order @f$q@f$ 
%> (see <code>phiTaylorPhy()</code> for properties of multi-indices).
%> We determine the correction factor @f$\alpha_{ke}^{(q)}@f$ for each order 
%> @f$q \le p@f$ (@f$p@f$ being the polynomial degree of the solution)
%> by computing correction factors using the linear limiter for all
%> linear reconstructions of derivatives of order @f$q-1@f$,
%> @f[
%>  \forall\,\mathbf{a} \in \mathcal{A}_{q-1}\,, \quad  
%>  c_{k,\mathbf{a},i} := C_{k,I(\mathbf{a})}\,\phi_{k1}(\mathbf{x}_{ki})
%>  + C_{k,I(\mathbf{a}+[1,0]^T)}\,\phi_{k2}(\mathbf{x}_{ki})
%>  + C_{k,I(\mathbf{a}+[0,1]^T)}\,\phi_{k3}(\mathbf{x}_{ki})
%>  \quad\text{on}~T_k\in\mathcal{T}_h\,,
%> @f]
%> where the indices of the corresponding degrees of freedom are given by 
%> @f$I(\mathbf{a})@f$ (cf. <code>mult2ind()</code>) and the first 
%> @f$x^1@f$- @f$I(\mathbf{a}+[1,0]^T)@f$ and @f$x^2@f$-derivatives 
%> @f$I(\mathbf{a}+[0,1]^T)@f$. 
%> The correction factor is defined as
%> @f[
%> \alpha_{ke}^{(q)} = \min_{\mathbf{a}\in\mathcal{A}_{q-1}} \, \alpha_{k\mathbf{a}}^{(q)} \,,
%> \qquad\text{with}\quad
%> \alpha_{k\mathbf{a}}^{(q)} \;:=\; \min_i \left\{
%> \begin{array}{llr}
%>  (c_{k,\mathbf{a},i}^\mathrm{max} - c_{k,\mathbf{a},c})\big/(c_{k,\mathbf{a},i} - c_{k,\mathbf{a},c})&\text{if}\;&c_{k,\mathbf{a},i} > c_{k,\mathbf{a},i}^\mathrm{max}\\
%>  1                  &\text{if}\;&c_{k,\mathbf{a},i}^\mathrm{min}\le c_{k,\mathbf{a},i}\le c_{k,\mathbf{a},i}^\mathrm{max}\\
%>  (c_{k,\mathbf{a},i}^\mathrm{min} - c_{k,\mathbf{a},c})\big/(c_{k,\mathbf{a},i} - c_{k,\mathbf{a},c})&\text{if}\;&c_{k,\mathbf{a},i} < c_{k,\mathbf{a},i}^\mathrm{min}
%> \end{array}\right\}\;,
%> @f]
%> where @f$c_{k,\mathbf{a},i}^\mathrm{min}@f$, @f$c_{k,\mathbf{a},i}^\mathrm{max}@f$
%> are minimum- and maximum values of above reconstruction as determined by
%> <code>computeMinMaxV0TElementPatch()</code>, @f$c_{k,\mathbf{a},i}@f$
%> are the function values in the vertices 
%> @f$\mathbf{x}_{ki}, i\in\{1,2,3\}@f$, and @f$c_{k,\mathbf{a},c}@f$ 
%> is the function value in the centroid @f$\mathbf{x}_{kc}@f$ of element
%> @f$T_k@f$.
%>
%> To avoid the loss of accuracy at smooth extrema, the lower order
%> derivatives should be limited by a factor not exceeding that of the 
%> higher order derivatives, since lower orders are typically smoother.
%> Beginning with the highest-order degrees of freedom, we compute the 
%> correction factors
%> @f[
%> \forall q \ge 1\,\quad  \alpha_{ke}^{(q)} \;\:=\; \max_{q \le d \le p} \alpha_{ke}^{(d)}  \,.
%> @f]
%>
%> The limited solution becomes then
%> @f[
%> c_h(\mathbf{x}) = C_{k,1} \,\phi_{k1} 
%> + \alpha_{ke}^{(1)} C_{k,2} \,\phi_{k2}(\mathbf{x})
%> + \alpha_{ke}^{(1)} C_{k,3} \,\phi_{k3}(\mathbf{x})
%> + \sum_{i=4}^{N} \alpha_{ke}^{(|\mathbf{a}_i|)} C_{ki}\, \phi_{ki}(\mathbf{x}) \,.
%> @f]
%>
%> For details to the hierarchical vertex based limiter, see:
%> 
%> D. Kuzmin, Slope limiting for discontinuous Galerkin approximations
%> with a possibly non-orthogonal Taylor basis, International Journal for
%> Numerical Methods in Fluids 71 (9) (2013) 1178–1190. doi:10.1002/fld.3707.
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  dataTaylor The representation matrix of the unlimited function 
%>                    @f$c_h@f$. @f$[K \times N]@f$
%> @param  markV0TbdrD <code>logical</code> arrays that mark each triangles
%>                    (Dirichlet boundary) vertices on which additional
%>                    function values should be considered during the 
%>                    slope limiting routine. @f$[K \times 3]@f$
%> @param  dataV0T    The function values for (Dirichlet boundary) vertices
%>                    specified by <code>markV0TbdrD</code>. @f$[K \times 3]@f$
%> @param  basesOnQuad  A struct containing precomputed values of (Taylor) basis
%>                      functions on quadrature points. Must provide at
%>                      least phiTaylorV0T.
%> @retval dataTaylorLim   The representation matrix of the limited function
%>                    @f$\mathsf{\Phi}^\mathrm{Taylor}c_h@f$. @f$[K \times N]@f$
%> @retval minMaxV0T   Two matrices with minimum and maximum centroid values,
%>                     respectively, of the patch of elements surrounding each
%>                     vertex of each element as computed by 
%>                     <code>computeMinMaxV0TElementPatch()</code>
%>                     @f$[2 \times 1 \mathrm{cell}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Florian Frank, Balthasar Reuter, Vadym Aizinger
%>
%> @author Hennes Hajduk, 2016
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
function [dataTaylorLim, minMaxV0T] = applySlopeLimiterTaylorHierarchicalVertex(g, dataTaylor, markV0TbdrD, dataV0T, basesOnQuad)
% Number of elements and polynomial degree
[K, N] = size(dataTaylor);
p = (sqrt(8*N+1)-3)/2;

% Check function arguments that are directly used
validateattributes(dataTaylor, {'numeric'}, {'size', [g.numT NaN]}, mfilename, 'dataTaylor');
assert(size(dataTaylor, 2) >= 3, 'Number of local degrees of freedom in dataTaylor does not correspond to p>=1')
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
validateattributes(basesOnQuad.phiTaylorV0T, {'numeric'}, {'size', [K 3 N]}, mfilename, 'phiTaylorV0T');

% Initialize limited coefficients and limiter of previous order
dataTaylorLim = zeros(size(dataTaylor));
dataTaylorLim(:, 1) = dataTaylor(:, 1);
alpha = zeros(K, 1);

% Determine limiter parameter for each order, beginning with highest order
for ord = p : -1 : 1
  alphaOrd = ones(K, 1);
  indDOF = ord*(ord+1)/2 + 1 : (ord+1)*(ord+2)/2;
  
  % Number of necessary reconstructions is equal to polynomial degree
  for i = 1 : ord
    % Determine multi-indices for linear reconstruction, shift them by
    %  multi-index of the centroid, and compute corresponding linear indices
    mult = bsxfun(@plus, [ord - i, i - 1], multiindex(1));
    ind = mult2ind(mult);
    
    % Compute limiter parameters for each element
    valV0T = computeFuncDiscAtPoints(dataTaylor(:, ind), basesOnQuad.phiTaylorV0T(:, :, 1:3));
    if ord > 1
      [alphaTmp, minMaxV0T] = computeVertexBasedLimiter(g, dataTaylor(:, ind(1)), valV0T, markV0TbdrD, valV0T);
    else
      [alphaTmp, minMaxV0T] = computeVertexBasedLimiter(g, dataTaylor(:, ind(1)), valV0T, markV0TbdrD, dataV0T);
    end
    alphaOrd = min(alphaOrd, alphaTmp);
  end %for
  
  % Determine limiter as hierarchical coarsening scheme
  alpha = max(alpha, alphaOrd);
  
  % Apply limiter
  dataTaylorLim(:, indDOF) = bsxfun(@times, alpha, dataTaylor(:, indDOF));
end %for
end % function

