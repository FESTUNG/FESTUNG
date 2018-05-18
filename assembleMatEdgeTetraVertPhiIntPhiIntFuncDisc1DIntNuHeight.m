% Assembles a matrix containing integrals over vertical edges of products of
% two two-dimensional basis function with a discrete 1D function, all from the
% interior of the element, and the x1-component of the edge normal, divided by 
% the smoothed mesh height.

%===============================================================================
%> @file
%>
%> @brief Assembles a matrix containing integrals over vertical edges of products
%>        of two two-dimensional basis function with a discrete 1D function, all
%>        from the interior of the element, and the @f$x^1@f$-component of the 
%>        edge normal, divided by the smoothed mesh height @f$H_\mathrm{s}@f$.
%===============================================================================
%>
%> @brief Assembles a matrix containing integrals over vertical edges of products
%>        of two two-dimensional basis function with a discrete 1D function, all
%>        from the interior of the element, and the @f$x^1@f$-component of the 
%>        edge normal, divided by the smoothed mesh height @f$H_\mathrm{s}@f$.
%>
%> The matrix @f$\check{\mathsf{P}}_\mathrm{bdr} \in 
%>      \mathbb{R}^{KN\times KN}@f$
%> is essentially the same as the matrix assembled by 
%> swe_2dv/assembleMatEdgeTetraVertPhiPhiFuncDisc1DNuHeight.m restricted to
%> diagonal blocks.
%> These are defined as
%> @f[
%> \left[\check{\mathsf{P}}_\mathrm{bdr}\right]_{(k-1)N+i,(k-1)N+j} :=
%>   \frac{1}{2} \sum_{E_{kn}\in\partial T_k \cap \mathcal{E}_\Omega^\mathrm{v}}
%>   \nu_{kn}^1 \sum_{l=1}^{\overline{N}} H_{\overline{k}l}
%>   \int_{E_{kn}} \frac{1}{H_\mathrm{s}} \varphi_{ki}\,\phi_{\overline{k}l}\,
%>   \varphi_{kj} \, \mathrm{d}\sigma \,,
%> @f]
%> with @f$\nu_{kn}^1@f$ the @f$x^1@f$-component (@f$m\in\{1,2\}@f$) of the edge
%> normal and @f$H_{\overline{k}l}@f$ the @f$(k,l)@f$-th entry of the representation
%> matrix of the discrete function @f$h_\Delta@f$.
%>
%> @param  g2D        The lists describing the geometric and topological 
%>                    properties of a quadrilateral triangulation (see 
%>                    <code>domainRectTrap()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  g1D        The lists describing the geometric and topological 
%>                    properties of a matching 1D triangulation (see 
%>                    <code>generateGridData1D()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param dataDisc1D  A representation matrix of the discrete function 
%>                    @f$h_\Delta(x^1)@f$, e.g., as computed by 
%>                    <code>projectFuncCont2DataDisc1D()</code>
%>                    @f$[\overline{K} \times \overline{N}]@f$
%> @param heightV0T1D The smoothed mesh height in each node of the 1D grid
%>                    @f$[\overline{K} \times 2]@f$
%> @param  markE0T    <code>logical</code> arrays that mark each elements
%>                    edges on which the matrix blocks should be
%>                    assembled @f$[K \times 4]@f$
%> @param refEdgePhiIntPhiIntPhi1DInt  Local matrix 
%>                    @f$\hat{\mathsf{Q}}^\text{diag}@f$ as provided
%>                    by <code>integrateRefEdgeTetraPhiIntPhiIntPhi1DInt()</code>.
%>                    @f$[N \times N \times N \times 4]@f$
%>
%> @retval ret        The assembled matrix @f$[KN \times KN]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
%>
%> @author Balthasar Reuter, 2017
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
function ret = assembleMatEdgeTetraVertPhiIntPhiIntFuncDisc1DIntNuHeight(g2D, g1D, dataDisc1D, heightV0T1D, markE0T, refEdgePhiIntPhiIntPhi1DInt)
K = g2D.numT;
[N, ~, barN, ~] = size(refEdgePhiIntPhiIntPhi1DInt);
ret = sparse(K*N, K*N);
dataDisc2D = g1D.markT2DT * dataDisc1D;
for n = 3 : 4
  areaNuE0Tint = markE0T(:,n) .* g2D.areaE0T(:,n) .* g2D.nuE0T(:,n,1);
  areaNuE0THeightint = areaNuE0Tint ./ (g1D.markT2DT * heightV0T1D(:,5-n));
  for l = 1 : barN
    ret = ret + kron(spdiags(areaNuE0THeightint .* dataDisc2D(:,l), 0, K, K), refEdgePhiIntPhiIntPhi1DInt(:,:,l,n));
  end % for l
end  % for n
end  % function