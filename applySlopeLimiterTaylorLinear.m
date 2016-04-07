% Applies the linear vertex based slope limiter to a discrete 
% function given in Taylor basis.
%
%===============================================================================
%> @file applySlopeLimiterTaylorLinear.m
%>
%> @brief Applies the linear vertex based slope limiter to a discrete 
%>        function given in Taylor basis.
%===============================================================================
%>
%> @brief Applies the linear vertex based slope limiter to a discrete 
%>        function given in Taylor basis.
%>
%> The goal is to determine the maximum admissible slope for a linear
%> reconstruction of the form
%> @f[
%> c_h(\mathbf{x}) = c_{k\mathrm{c}} + 
%>   \alpha_{ke} \, \nabla c (\mathbf{x}_{k\mathrm{c}}) \cdot 
%>   (\mathbf{x} - \mathbf{x}_{k\mathrm{c}}) \,,
%> \qquad 0\le \alpha_e\le 1\,,\quad \mathbf{x} \in T_k\,,
%> @f]
%> where we abbreviated the function value 
%> @f$c_{k\mathrm{c}} := c_h(\mathbf{x}_{k\mathrm{c}})@f$ in the centroid
%> @f$\mathbf{x}_{k\mathrm{c}}@f$.
%> The correction factor @f$\alpha_{ke}@f$ is chosen such that above
%> reconstruction is bounded in all vertices @f$\mathbf{x}_{ki} \in T_k@f$
%> by the minimum and maximum centroid values of all elements containing
%> @f$\mathbf{x}_{ki}@f$, that is
%> @f[
%> \forall T_k \in \mathcal{T}_h \,, \forall i\in\{1,2,3\} \,,\quad 
%> c_{ki}^\mathrm{min} \le c_h(\mathbf{x}_{ki}) \le c_{ki}^\mathrm{max}\,.
%> @f]
%> See <code>computeMinMaxV0TElementPatch</code>
%> for details to the bounds.
%>
%> The correction factor @f$\alpha_{ke}@f$ is defined as
%> @f[
%> \forall T_k\in\mathcal{T}_h\,,\quad  
%> \alpha_{ke} := \min_{i\,\in\, \{1, 2, 3\}} \left\{
%> \begin{array}{llr}
%>  (c_{ki}^\mathrm{max} - c_{k\mathrm{c}})\big/(c_{ki} - c_{k\mathrm{c}})&\text{if}\;&c_{ki} > c_{ki}^\mathrm{max}\\
%>  1                  &\text{if}\;&c_{ki}^\mathrm{min}\le c_{ki}\le c_{ki}^\mathrm{max}\\
%>  (c_{ki}^\mathrm{min} - c_{k\mathrm{c}})\big/(c_{ki} - c_{k\mathrm{c}})&\text{if}\;&c_{ki} < c_{ki}^\mathrm{min}
%> \end{array}\right\}\;,
%> @f]
%> where @f$c_{ki} := c_{k\mathrm{c}} + \nabla c (\mathbf{x}_{k\mathrm{c}}) 
%> \cdot (\mathbf{x}_{ki} - \mathbf{x}_{k\mathrm{c}})@f$ 
%> is the unconstrained linear reconstruction in @f$\mathbf{x}_{ki}@f$.
%>
%> The limited solution becomes then
%> @f[
%> c_h(\mathbf{x}) = C_{k,1} \,\phi_{k1} 
%> + \alpha_{ke} ( C_{k,2} \,\phi_{k2}(\mathbf{x}) + C_{k,3} \,\phi_{k3}(\mathbf{x}) )
%> + \delta_{\alpha_{ke} = 1} \sum_{i=4}^{N} C_{ki}\, \phi_{ki}(\mathbf{x}) \,.
%> @f]
%>
%> For details to the linear vertex based limiter, see:
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
%> @retval dataTaylorLim   The representation matrix of the limited function
%>                    @f$\mathsf{\Phi}^\mathrm{Taylor}c_h@f$. @f$[K \times N]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function dataTaylorLim = applySlopeLimiterTaylorLinear(g, dataTaylor, markV0TbdrD, dataV0T)
global gPhiTaylorV0T

% Check function arguments that are directly used
validateattributes(dataTaylor, {'numeric'}, {'size', [g.numT NaN]}, mfilename, 'dataTaylor');
assert(size(dataTaylor, 2) >= 3, 'Number of local degrees of freedom in dataTaylor does not correspond to p>=1')
validateattributes(gPhiTaylorV0T, {'numeric'}, {'size', [K 3 N]}, mfilename, 'gPhiTaylorV0T');

%% Limit first order derivative terms
% Compute limiter parameter for each vertex
alphaE = computeVertexBasedLimiter(g, dataTaylor(:, 1), computeFuncDiscAtPoints(dataTaylor(:, 1:3), gPhiTaylorV0T(:,:,1:3)), markV0TbdrD, dataV0T);

% Apply limiter to first order terms, set all higher order terms to zero in
% the case of limiting
dataTaylorLim = [dataTaylor(:,1), bsxfun(@times, alphaE, dataTaylor(:, 2:3)), bsxfun(@times, alphaE == 1, dataTaylor(:, 4:end))];
end % function

