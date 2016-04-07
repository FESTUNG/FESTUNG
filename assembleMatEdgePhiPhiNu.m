% Assembles two matrices containing integrals over interior edges of products of
% two basis functions with a component of the edge normal.
%
%===============================================================================
%> @file assembleMatEdgePhiPhiNu.m
%>
%> @brief Assembles two matrices containing integrals over interior edges of 
%>        products of two basis functions with a component of the edge normal.
%===============================================================================
%>
%> @brief Assembles matrices @f$\mathsf{{Q}}^m, m \in \{1,2\}@f$
%>        containing integrals over edges of products of two basis functions 
%>        with a component of the edge normal.
%>
%> The matrices @f$\mathsf{{Q}}^m \in \mathbb{R}^{KN\times KN}@f$
%> consist of two kinds of contributions: diagonal blocks and off-diagonal 
%> blocks. Diagonal blocks are defined as 
%> @f[
%> [\mathsf{{Q}}^m]_{(k-1)N+i,(k-1)N+j} = \frac{1}{2}
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_\Omega}
%>  \nu_{kn}^m  \int_{E_{kn}} \varphi_{ki} \varphi_{kj} \mathrm{d}s \,.
%> @f]
%> and off-diagonal blocks are defined as
%> @f[
%> [\mathsf{{Q}}^m]_{(k^--1)N+i,(k^+-1)N+j} = \frac{1}{2}
%>  \nu_{k^-n^-}^m \int_{E_{k^-n^-}} 
%>   \varphi_{k^-i} \varphi_{k^+j} \mathrm{d}s \,.
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
%>     \nu^m_{1n} | E_{1n} | &   & \\
%>     & ~\ddots~ & \\
%>     &          & \nu^m_{Kn} | E_{Kn} |
%>   \end{bmatrix} \otimes [\hat{\mathsf{Q}}^\mathrm{diag}]_{:,:,n}\;,
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
%>     \nu_{1n^-}^m |E_{1n^-}| & \dots & \nu_{1n^-}^m |E_{1n^-}| \\
%>     \vdots & & \vdots \\
%>     \nu_{Kn^-}^m |E_{Kn^-}| & \dots & \nu_{Kn^-}^m |E_{Kn^-}| 
%>   \end{bmatrix} \otimes 
%> [\hat{\mathsf{Q}}^\mathrm{offdiag}]_{:,:,n^-,n^+}\;,
%> @f]
%> where @f$\delta_{E_{kn}\in\mathcal{E}_\Omega}, \delta_{E_{in^-}=E_{jn^+}}@f$ 
%> denote the Kronecker delta, @f$\circ@f$ denotes the Hadamard product, and 
%> @f$\otimes@f$ denotes the Kronecker product.
%>
%> The entries of matrix 
%> @f$\hat{\mathsf{Q}}^\mathrm{diag}\in\mathbb{R}^{N\times N\times3}@f$
%> are given by
%> @f[
%> [\hat{\mathsf{Q}}^\mathrm{diag}]_{i,j,n} =
%>   \int_0^1 \hat{\varphi}_i \circ \hat{\mathbf{\gamma}}_n(s) 
%>   \hat{\varphi}_j\circ \hat{\mathbf{\gamma}}_n(s) \mathrm{d}s \,,
%> @f]
%> where the mapping @f$\hat{\mathbf{\gamma}}_n@f$ is defined in 
%> <code>gammaMap()</code>. The entries of matrix
%> @f$\hat{\mathsf{Q}}^\mathrm{offdiag} \in 
%>    \mathbb{R}^{N\times N\times 3\times 3}@f$ are defined as
%> @f[
%> [\hat{\mathsf{Q}}^\mathrm{offdiag}]_{i,j,n^-,n^+} =
%>   \int_0^1 \hat{\varphi}_i \circ \hat{\mathbf{\gamma}}_{n^-}(s) 
%>   \hat{\varphi}_j\circ \hat{\mathbf{\vartheta}}_{n^-n^+} \circ
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
%> @param refEdgePhiIntPhiInt  Local matrix 
%>                    @f$\hat{\mathsf{Q}}^\text{diag}@f$ as provided
%>                    by <code>integrateRefEdgePhiIntPhiInt()</code>.
%>                    @f$[N \times N \times 3]@f$
%> @param refEdgePhiIntPhiExt Local matrix 
%>                    @f$\hat{\mathsf{Q}}^\text{offdiag}@f$ as provided
%>                    by <code>integrateRefEdgePhiIntPhiExt()</code>.
%>                    @f$[N \times N \times 3 \times 3]@f$
%> @param areaNuE0Tint (optional) argument to provide precomputed values
%>                    for the products of <code>markE0Tint</code>,
%>                    <code>g.areaE0T</code>, and <code>g.nuE0T</code>
%>                    @f$[3 \times 2 \text{ cell}]@f$
%> @retval ret        The assembled matrices @f$[2 \times 1 \text{ cell}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2015 Florian Frank, Balthasar Reuter, Vadym Aizinger
%> Modified by Hennes Hajduk, 2016-04-06
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
function ret = assembleMatEdgePhiPhiNu(g, markE0Tint, refEdgePhiIntPhiInt, refEdgePhiIntPhiExt, areaNuE0Tint)
K = g.numT;  N = size(refEdgePhiIntPhiInt, 1);

% Check function arguments that are directly used
validateattributes(markE0Tint, {'logical'}, {'size', [K 3]}, mfilename, 'markE0Tint');
validateattributes(refEdgePhiIntPhiInt, {'numeric'}, {'size', [N N 3]}, mfilename, 'refEdgePhiIntPhiInt');
validateattributes(refEdgePhiIntPhiExt, {'numeric'}, {'size', [N N 3 3]}, mfilename, 'refEdgePhiIntPhiExt');

% Assemble matrices
ret = cell(2, 1); ret{1} = sparse(K*N, K*N);  ret{2} = sparse(K*N, K*N);
for nn = 1 : 3
  if isfield(g, 'areaNuE0TE0T')
    for np = 1 : 3
      ret{1} = ret{1} + 0.5 * kron(g.areaNuE0TE0T{nn,np,1}, refEdgePhiIntPhiExt(:,:,nn,np));
      ret{2} = ret{2} + 0.5 * kron(g.areaNuE0TE0T{nn,np,2}, refEdgePhiIntPhiExt(:,:,nn,np));
    end % for
  else
    Qkn = 0.5 * g.areaE0T(:,nn);
    for np = 1 : 3
      ret{1} = ret{1} + kron(bsxfun(@times, g.markE0TE0T{nn,np}, Qkn .* g.nuE0T(:,nn,1)), refEdgePhiIntPhiExt(:,:,nn,np));
      ret{2} = ret{2} + kron(bsxfun(@times, g.markE0TE0T{nn,np}, Qkn .* g.nuE0T(:,nn,2)), refEdgePhiIntPhiExt(:,:,nn,np));
    end % for
  end % if
  if nargin > 4
    ret{1} = ret{1} + kron(spdiags(0.5 * areaNuE0Tint{nn,1}, 0,K,K), refEdgePhiIntPhiInt(:,:,nn));
    ret{2} = ret{2} + kron(spdiags(0.5 * areaNuE0Tint{nn,2}, 0,K,K), refEdgePhiIntPhiInt(:,:,nn));
  else
    if isfield(g, 'areaNuE0T')
      ret{1} = ret{1} + kron(spdiags(0.5 * g.areaNuE0T{nn,1} .* markE0Tint(:, nn), 0,K,K), refEdgePhiIntPhiInt(:,:,nn));
      ret{2} = ret{2} + kron(spdiags(0.5 * g.areaNuE0T{nn,2} .* markE0Tint(:, nn), 0,K,K), refEdgePhiIntPhiInt(:,:,nn));
    else
      if ~isfield(g, 'areaNuE0TE0T')
        Qkn = Qkn .* markE0Tint(:,nn);
      else
        Qkn = 0.5 * g.areaE0T(:,nn) .* markE0Tint(:, nn);
      end % if
      ret{1} = ret{1} + kron(spdiags(Qkn .* g.nuE0T(:,nn,1), 0,K,K), refEdgePhiIntPhiInt(:,:,nn));
      ret{2} = ret{2} + kron(spdiags(Qkn .* g.nuE0T(:,nn,2), 0,K,K), refEdgePhiIntPhiInt(:,:,nn));
    end
  end % if
end % for
end % function
