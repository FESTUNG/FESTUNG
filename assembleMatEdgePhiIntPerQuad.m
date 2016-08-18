% Assembles three matrices containing evaluations of one basis function in 
% quadrature points for each of the three local edges of each element multiplied
% with the corresponding quadrature weight.

%===============================================================================
%> @file assembleMatEdgePhiIntPerQuad.m
%>
%> @brief Assembles three matrices containing evaluations of one basis function 
%>        in quadrature points for each of the three local edges of each element
%>        multiplied with the corresponding quadrature weight.
%===============================================================================
%>
%> @brief Assembles three matrices containing evaluations of one basis function 
%>        in quadrature points for each of the three local edges of each element
%>        multiplied with the corresponding quadrature weight.
%>
%> The matrix @f$\mathsf{{Q}_n \in \mathbb{R}^{KN\times KR}@f$ (R is the 
%> number of quadrature points and weights.) is block diagonal and defined as 
%> @f[
%> [\mathsf{{Q}}_n]_{(k-1)N+i,(k-1)R+j} = \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_N}
%>  \varphi_{ki}(q^r_{kn}) w^r_{kn} \,.
%> @f]
%> with q^r_{kn}, w^r_{kn} the quadrature points and weights of edge n of element k.
%> All other entries are zero.
%> To allow for vectorization, the assembly is reformulated as
%> @f[
%> \mathsf{{Q}}^m_n = 
%>   \begin{bmatrix}
%>     \delta_{E_{1n}\in\mathcal{E}_\mathrm{N}} &   & \\
%>     & ~\ddots~ & \\
%>     &          & \delta_{E_{Kn}\in\mathcal{E}_\mathrm{N}}
%>   \end{bmatrix} \circ \begin{bmatrix}
%>     | E_{1n} | &   & \\
%>     & ~\ddots~ & \\
%>     &          & | E_{Kn} |
%>   \end{bmatrix} 
%>  \otimes [\hat{\mathsf{{S}}}]_{:,:,n}\;,
%> @f]
%> where @f$\delta_{E_{kn}\in\mathcal{E}_\mathrm{N}}@f$ denotes the Kronecker 
%> delta, @f$\circ@f$ denotes the Hadamard product, and @f$\otimes@f$ denotes 
%> the Kronecker product.
%>
%> The entries of matrix 
%> @f$\hat{\mathsf{{S}}}\in\mathbb{Q}^{N\times R \times 3}@f$
%> are given by
%> @f[
%> [\hat{\mathsf{{S}}}]_{i,r,n} =
%>   \hat{\varphi}_i \circ \hat{\mathbf{\gamma}}_n(q_r) w_r\,,
%> @f]
%> where the mapping @f$\hat{\mathbf{\gamma}}_n@f$ is defined in 
%> <code>gammaMap()</code>.
%>
%> It is essentially the same as the diagonal part of
%> <code>assembleMatEdgePhiNuPerQuad()</code>.
%>
%> @param g           The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param markE0Tbdr  <code>logical</code> arrays that mark each triangles
%>                    (boundary) edges on which the matrix blocks should be
%>                    assembled @f$[K \times 3]@f$
%> @param refEdgePhiIntPerQuad  Local matrix 
%>                    @f$\hat{\mathsf{S}}@f$ as provided
%>                    by <code>integrateRefEdgePhiIntPerQuad()</code>.
%>                    @f$[N \times R \times  3]@f$
%> @param areaE0Tbdr (optional) argument to provide precomputed values
%>                    for the products of <code>markE0Tbdr</code> and
%>                    <code>g.areaE0T</code>,
%>                    @f$[3 \times 1 \text{ cell}]@f$
%> @retval ret        The assembled matrices @f$[3 \times 1 \text{ cell}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Hennes Hajduk, Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function ret = assembleMatEdgePhiIntPerQuad(g, markE0Tbdr, refEdgePhiIntPerQuad, areaE0Tbdr)
K = g.numT;
ret = cell(3,1);
for n = 1 : 3
  if nargin > 3
    ret{n} = 0.5 * kron( spdiags(areaE0Tbdr{n}, 0, K, K), refEdgePhiIntPerQuad(:,:,n) );
  else
    ret{n} = 0.5 * kron( spdiags(markE0Tbdr(:,n) .* g.areaE0T(:,n), 0, K, K), refEdgePhiIntPerQuad(:,:,n) );
  end % if
end % for
end % function
