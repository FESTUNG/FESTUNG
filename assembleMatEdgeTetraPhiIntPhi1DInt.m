% Assembles a matrix containing integrals over edges of products of
% a two-dimensional and a one-dimensional basis function.

%===============================================================================
%> @file
%>
%> @brief Assembles a matrix containing integrals over edges of 
%>        products of a two-dimensional and a one-dimensional basis function.
%===============================================================================
%>
%> @brief Assembles a matrix containing integrals over edges of 
%>        products of a two-dimensional and a one-dimensional basis function.
%>
%> The matrix @f$\mathsf{{Q}}^{\mathrm{diag}} \in 
%>      \mathbb{R}^{KN\times \overline{K}\overline{N}}@f$ is defined as 
%> @f[
%> [\mathsf{{Q}}^{\mathrm{diag}}]_{(k-1)N+i,(\overline{k}-1)\overline{N}+j} = 
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_\Omega}
%>   \int_{E_{kn}} \varphi_{ki} \phi_{\overline{k}j} \mathrm{d}s \,.
%> @f]
%> To allow for vectorization, the assembly is reformulated as
%> @f[
%> \mathsf{Q}^{\mathrm{diag}} = \sum_{n=1}^4
%>   \begin{bmatrix}
%>     \delta_{E_{1n}\in\mathcal{E}_\Omega} &   & \\
%>     & ~\ddots~ & \\
%>     &          & \delta_{E_{Kn}\in\mathcal{E}_\Omega}
%>   \end{bmatrix} \circ \begin{bmatrix}
%>     | E_{1n} | &   & \\
%>     & ~\ddots~ & \\
%>     &          & | E_{Kn} |
%>   \end{bmatrix} \otimes [\hat{\mathsf{Q}}^\mathrm{diag}]_{:,:,n}\;,
%> @f]
%> where @f$\delta_{E_{kn}\in\mathcal{E}_\Omega}@f$ denotes the Kronecker delta, 
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
%> @param coefE0T     (optional) Coefficient to be applied to each edge
%>                    integral. Defaults to <code>g2D.areaE0T</code>
%>                    @f$[K\times 4]@f$
%> @retval ret        The assembled matrix @f$[2 \times 1 \text{ cell}]@f$
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
function ret = assembleMatEdgeTetraPhiIntPhi1DInt(g2D, g1D, markE0T, refEdgePhiIntPhi1DInt, coefE0T)
if nargin < 5
  coefE0T = g2D.areaE0T;
end % if
K = g2D.numT; barK = g1D.numT;
[N, barN, ~] = size(refEdgePhiIntPhi1DInt);
ret = sparse(K*N, barK*barN);
for n = 1 : 4
  ret = ret + kron(bsxfun(@times, g1D.markT2DT, markE0T(:,n) .* coefE0T(:,n)), refEdgePhiIntPhi1DInt(:,:,n));
end  % for n
end  % function