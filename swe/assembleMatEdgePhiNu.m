% Assembles two matrices containing integrals over interior edges of products of
% a basis function with a component of the edge normal divided by the element 
% area.

%===============================================================================
%> @file
%>
%> @brief Assembles two matrices containing integrals over interior edges of 
%>        products ofa basis function with a component of the edge normal 
%>        divided by the element area.
%===============================================================================
%>
%> @brief Assembles two matrices @f$\mathsf{{Q}}^m, m \in \{1,2\}@f$containing 
%>        integrals over interior edges of  products ofa basis function with a 
%>        component of the edge normal divided by the element area.
%>
%> The matrices @f$\mathsf{{Q}}^m \in \mathbb{R}^{K\times KN}@f$
%> consist of two kinds of contributions: diagonal blocks and off-diagonal 
%> blocks. Diagonal blocks are defined as 
%> @f[
%> [\mathsf{{Q}}^m]_{k,(k-1)N+j} = \frac{1}{2}
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_\Omega}
%>  \nu_{kn}^m / |T_k| \int_{E_{kn}} \varphi_{kj} \mathrm{d}s \,.
%> @f]
%> and off-diagonal blocks are defined as
%> @f[
%> [\mathsf{{Q}}^m]_{k^-,(k^+-1)N+j} = \frac{1}{2}
%>  \nu_{k^-n^-}^m / |T_k| \int_{E_{k^-n^-}} 
%>    \varphi_{k^+j} \mathrm{d}s \,.
%> @f]
%> with @f$\nu_{kn}@f$ the @f$m@f$-th component (@f$m\in\{1,2\}@f$) of the edge
%> normal.
%> Entries in off-diagonal blocks are only non-zero for pairs of triangles
%> @f$T_{k^-}, T_{k^+}@f$ with @f$\partial T_{k^-} \cap T_{k^+} \ne\emptyset@f$.
%> Note that the local edge index @f$n^-@f$ is given implicitly, since 
%> @f$\partial T_{k^-} \cap T_{k^+} \ne\emptyset@f$ consist of exactly one
%> edge @f$E_{k^-n^-} = E_{k^+n^+}@f$.
%> To allow for vectorization, the assembly is reformulated as
%> @f$\mathsf{Q}^m = \mathsf{Q}^{m,\mathrm{diag}} + 
%>    \mathsf{Q}^{m,\mathrm{offdiag}}@f$ with the blocks defined as
%> @f[
%> \mathsf{{Q}}^{m,\mathrm{diag}} = \frac{1}{2} \sum_{n=1}^3
%>   \begin{bmatrix}
%>     \delta_{E_{1n}\in\mathcal{E}_\Omega} &   & \\
%>     & ~\ddots~ & \\
%>     &          & \delta_{E_{Kn}\in\mathcal{E}_\Omega}
%>   \end{bmatrix} \circ \begin{bmatrix}
%>     \nu^m_{1n} | E_{1n} | / |T_1| &   & \\
%>     & ~\ddots~ & \\
%>     &          & \nu^m_{Kn} | E_{Kn} | / |T_K|
%>   \end{bmatrix} \otimes [\hat{\mathsf{Q}}^\mathrm{diag}]_{:,n}\;,
%> @f]
%> and
%> @f[
%> \mathsf{Q}^{m,\mathrm{offdiag}} = \frac{1}{2}\sum_{n^-=1}^3\sum_{n^+=1}^3
%>   \begin{bmatrix}
%>     0&\delta_{E_{1n^-} = E_{2n^+}}&\dots&\dots&\delta_{E_{1n^-}=E_{Kn^+}} \\
%>     \delta_{E_{2n^-} = E_{1n^+}}&0&\ddots& &\vdots \\
%>     \vdots & \ddots & \ddots & \ddots & \vdots \\
%>     \vdots & & \ddots & 0 & \delta_{E_{(K-1)n^-}=E_{Kn^+}} \\
%>     \delta_{E_{Kn^-} = E_{1n^+}}&\dots&\dots&\delta_{E_{Kn^-} = E_{(K-1)n^+}} &0
%>   \end{bmatrix} \circ \begin{bmatrix}
%>     \nu_{1n^-}^m |E_{1n^-}| / |T_1| & \dots & \nu_{1n^-}^m |E_{1n^-}| / |T_1| \\
%>     \vdots & & \vdots \\
%>     \nu_{Kn^-}^m |E_{Kn^-}| / |T_K| & \dots & \nu_{Kn^-}^m |E_{Kn^-}| / |T_K| 
%>   \end{bmatrix} \otimes 
%> [\hat{\mathsf{Q}}^\mathrm{offdiag}]_{:,n^-,n^+}\;,
%> @f]
%> where @f$\delta_{E_{kn}\in\mathcal{E}_\Omega}, \delta_{E_{in^-}=E_{jn^+}}@f$ 
%> denote the Kronecker delta, @f$\circ@f$ denotes the Hadamard product, and 
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
%> <code>gammaMap()</code>. The entries of matrix
%> @f$\hat{\mathsf{Q}}^\mathrm{offdiag} \in 
%>    \mathbb{R}^{N\times 3\times 3}@f$ are defined as
%> @f[
%> [\hat{\mathsf{Q}}^\mathrm{offdiag}]_{j,n^-,n^+} =
%>   \int_0^1 \hat{\varphi}_j\circ \hat{\mathbf{\vartheta}}_{n^-n^+} \circ
%>   \hat{\mathbf{\gamma}}_{n^-}(s) \mathrm{d}s \,,
%> @f]
%> with the mapping @f$\hat{\mathbf{\vartheta}}_{n^-n^+}@f$ as described in
%> <code>theta()</code>.
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
%> @param refEdgePhiExt Local matrix 
%>                    @f$\hat{\mathsf{Q}}^\text{offdiag}@f$ as provided
%>                    by <code>integrateRefEdgePhiExt()</code>.
%>                    @f$[N \times 3 \times 3]@f$
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
function ret = assembleMatEdgePhiNu(g, markE0Tint, refEdgePhiInt, refEdgePhiExt)
K = g.numT;  N = size(refEdgePhiInt, 1);

% Check function arguments that are directly used
validateattributes(markE0Tint, {'logical'}, {'size', [K 3]}, mfilename, 'markE0Tint');
validateattributes(refEdgePhiInt, {'numeric'}, {'size', [N 3]}, mfilename, 'refEdgePhiInt');
validateattributes(refEdgePhiExt, {'numeric'}, {'size', [N 3 3]}, mfilename, 'refEdgePhiExt');

% Assemble matrices
ret = cell(2, 1); 
ret{1} = sparse(K, K*N);  
ret{2} = sparse(K, K*N);
for nn = 1 : 3
  Qkn = 0.5 * g.areaE0T(:,nn) ./ g.areaT;
  for np = 1 : 3
    ret{1} = ret{1} + kron(bsxfun(@times, g.markE0TE0T{nn,np}, Qkn .* g.nuE0T(:,nn,1)), refEdgePhiExt(:,nn,np).');
    ret{2} = ret{2} + kron(bsxfun(@times, g.markE0TE0T{nn,np}, Qkn .* g.nuE0T(:,nn,2)), refEdgePhiExt(:,nn,np).');
  end % for
  Qkn = Qkn .* markE0Tint(:,nn);
  ret{1} = ret{1} + kron(spdiags(Qkn .* g.nuE0T(:,nn,1), 0,K,K), refEdgePhiInt(:,nn).');
  ret{2} = ret{2} + kron(spdiags(Qkn .* g.nuE0T(:,nn,2), 0,K,K), refEdgePhiInt(:,nn).');
end % for
end % function
