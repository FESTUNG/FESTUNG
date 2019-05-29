% Assembles two matrices containing integrals over boundary edges of products of
% a basis function with a component of the edge normal divided by the element 
% area.

%===============================================================================
%> @file
%>
%> @brief Assembles two matrices containing integrals over boundary edges of 
%>        products of a basis function with a component of the edge normal 
%>        divided by the element area.
%===============================================================================
%>
%> @brief Assembles two matrices @f$\mathsf{{Q}}^m_\mathrm{N}, m \in \{1,2\}@f$
%>        containing integrals over boundary edges of products of a basis 
%>        function with a component of the edge normal divided by the element 
%>        area.
%>
%> The matrices @f$\mathsf{{Q}}^m_\mathrm{N} \in \mathbb{R}^{K\times KN}@f$
%> are defined as 
%> @f[
%> [\mathsf{{Q}}^m_\mathrm{N}]_{k,(k-1)N+j} = 
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_\mathrm{N}}
%>  \nu_{kn}^m / |T_k| \int_{E_{kn}} \varphi_{kj} \mathrm{d}s \,.
%> @f]
%> with @f$\nu_{kn}@f$ the @f$m@f$-th component (@f$m\in\{1,2\}@f$) of the edge
%> normal.
%> To allow for vectorization, the assembly is reformulated as
%> \mathsf{{Q}}^m_\mathrm{N} = \sum_{n=1}^3
%>   \begin{bmatrix}
%>     \delta_{E_{1n}\in\mathcal{E}_\mathrm{N}} &   & \\
%>     & ~\ddots~ & \\
%>     &          & \delta_{E_{Kn}\in\mathcal{E}_\mathrm{N}}
%>   \end{bmatrix} \circ \begin{bmatrix}
%>     \nu^m_{1n} | E_{1n} | / |T_1| &   & \\
%>     & ~\ddots~ & \\
%>     &          & \nu^m_{Kn} | E_{Kn} | / |T_K|
%>   \end{bmatrix} \otimes [\hat{\mathsf{Q}}^\mathrm{diag}]_{:,n}\;,
%> @f]
%> @f]
%> where @f$\delta_{E_{kn}\in\mathcal{E}_\mathrm{N}}
%> denotes the Kronecker delta, @f$\circ@f$ denotes the Hadamard product, and 
%> @f$\otimes@f$ denotes the Kronecker product.
%>
%> The entries of matrix 
%> @f$\hat{\mathsf{Q}}^\mathrm{diag}\in\mathbb{R}^{N\times 3}@f$
%> are given by
%> @f[
%> [\hat{\mathsf{Q}}^\mathrm{diag}]_{j,n} =
%>   \int_0^1 \hat{\varphi}_j\circ \hat{\mathbf{\gamma}}_n(s) \mathrm{d}s \,,
%> @f]
%> where the mapping @f$\hat{\mathbf{\gamma}}_n@f$ is defined in 
%> <code>gammaMap()</code>. 
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0Tint <code>logical</code> arrays that mark each triangles
%>                    (interior) edges on which the matrix blocks should be
%>                    assembled @f$[K \times 3]@f$
%> @param refEdgePhiInt  Local matrix 
%>                    @f$\hat{\mathsf{Q}}^\text{diag}@f$ as provided
%>                    by <code>integrateRefEdgePhiInt()</code>.
%>                    @f$[N \times 3]@f$
%> @retval ret        The assembled matrices @f$[2 \times 1 \text{ cell}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2018 Florian Frank, Balthasar Reuter, Vadym Aizinger
%> author Hennes Hajduk, 2018
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
function ret = assembleMatEdgePhiIntNu(g, markE0Tbdr, refEdgePhiInt)
K = g.numT;  N = size(refEdgePhiInt, 1);

% Check function arguments that are directly used
validateattributes(markE0Tbdr, {'logical'}, {'size', [K 3]}, mfilename, 'markE0Tint');
validateattributes(refEdgePhiInt, {'numeric'}, {'size', [N 3]}, mfilename, 'refEdgePhiInt');

% Assemble matrices
ret = cell(2, 1); 
ret{1} = sparse(K, K*N);  
ret{2} = sparse(K, K*N);
for nn = 1 : 3
  Qkn = g.areaE0T(:,nn) .* markE0Tbdr(:,nn) ./ g.areaT;
  ret{1} = ret{1} + kron(spdiags(Qkn .* g.nuE0T(:,nn,1), 0,K,K), refEdgePhiInt(:,nn).');
  ret{2} = ret{2} + kron(spdiags(Qkn .* g.nuE0T(:,nn,2), 0,K,K), refEdgePhiInt(:,nn).');
end % for
end % function
