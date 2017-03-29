% Assembles two matrices containing integrals over edges of products of
% a two-dimensional and a one-dimensional basis function with a component of the
% edge normal.

%===============================================================================
%> @file assembleMatEdgeTetraPhiIntPhi1DIntNu.m
%>
%> @brief Assembles two matrices containing integrals over edges of 
%>        products of a two-dimensional and a one-dimensional basis function 
%>        with a component of the edge normal.
%===============================================================================
%>
%> @brief Assembles two matrices containing integrals over edges of 
%>        products of a two-dimensional and a one-dimensional basis function 
%>        with a component of the edge normal.
%>
%> The matrices @f$\mathsf{{Q}}^{m,\mathrm{diag}} \in 
%>      \mathbb{R}^{KN\times \overline{K}\overline{N}}@f$ are defined as 
%> @f[
%> [\mathsf{{Q}}^{m,\mathrm{diag}}]_{(k-1)N+i,(\overline{k}-1)\overline{N}+j} = 
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_\Omega}
%>  \nu_{kn}^m  \int_{E_{kn}} \varphi_{ki} \phi_{\overline{k}j} \mathrm{d}s \,.
%> @f]
%> with @f$\nu_{kn}@f$ the @f$m@f$-th component (@f$m\in\{1,2\}@f$) of the edge
%> normal.
%> To allow for vectorization, the assembly is reformulated as
%> \mathsf{{Q}}^{m,\mathrm{diag}} = \sum_{n=1}^4
%>   \begin{bmatrix}
%>     \delta_{E_{1n}\in\mathcal{E}_\Omega} &   & \\
%>     & ~\ddots~ & \\
%>     &          & \delta_{E_{Kn}\in\mathcal{E}_\Omega}
%>   \end{bmatrix} \circ \begin{bmatrix}
%>     \nu^m_{1n} | E_{1n} | &   & \\
%>     & ~\ddots~ & \\
%>     &          & \nu^m_{Kn} | E_{Kn} |
%>   \end{bmatrix} \otimes [\hat{\mathsf{Q}}^\mathrm{diag}]_{:,:,n}\;,
%> @f]
%> where @f$\delta_{E_{kn}\in\mathcal{E}_\Omega} denotes the Kronecker delta, 
%> @f$\circ@f$ denotes the Hadamard product, and 
%> @f$\otimes@f$ denotes the Kronecker product.
%>
%> The entries of matrix 
%> @f$\hat{\mathsf{Q}}^\mathrm{diag}\in\mathbb{R}^{N\times \overline{N}\times4}@f$
%> are given by
%> @f[
%> [\hat{\mathsf{Q}}^\mathrm{diag}]_{i,j,n} =
%>   \int_0^1 \hat{\varphi}_i \circ \hat{\mathbf{\gamma}}_n(s) 
%>   \hat{\phi}_j\circ [\hat{\mathbf{\gamma}}_n(s)]_1 \mathrm{d}s \,,
%> @f]
%> where the mapping @f$\hat{\mathbf{\gamma}}_n@f$ is defined in 
%> <code>gammaMapTetra()</code>.
%>
%> @param  g2D        The lists describing the geometric and topological 
%>                    properties of a quadrilateral triangulation (see 
%>                    <code>domainRectTrap()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  g1D        The lists describing the geometric and topological 
%>                    properties of a matching 1D triangulation (see 
%>                    <code>generateGridData1D()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0T    <code>logical</code> arrays that mark each triangles
%>                    edges on which the matrix blocks should be
%>                    assembled @f$[K \times 3]@f$
%> @param refEdgePhiIntPhi1DInt  Local matrix 
%>                    @f$\hat{\mathsf{Q}}^\text{diag}@f$ as provided
%>                    by <code>integrateRefEdgeTetraPhiIntPhi1DInt()</code>.
%>                    @f$[N \times N \times 4]@f$
%> @retval ret        The assembled matrices @f$[2 \times 1 \text{ cell}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
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
function ret = assembleMatEdgeTrapPhiIntPhi1DIntNu(g2D, g1D, markE0T, refEdgePhiIntPhi1DInt)
K = g2D.numT; barK = g1D.numT;
[N, barN, ~] = size(refEdgePhiIntPhi1DInt);
ret = { sparse(K*N, barK*barN), sparse(K*N, barK*barN) };
for n = 1 : 4
  areaE0Tint = markE0T(:,n) .* g2D.areaE0T(:,n);
  for m = 1 : 2
    areaNuE0Tint = areaE0Tint .* g2D.nuE0T(:,n,m);
    ret{m} = ret{m} + kron(bsxfun(@times, g1D.markT2DT, areaNuE0Tint), refEdgePhiIntPhi1DInt(:,:,n));
  end % for m
end  % for n
end  % function