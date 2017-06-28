% Assembles a matrix containing integrals over interior edges of products of two 
% basis functions.

%===============================================================================
%> @file assembleMatEdgeTetraPhiPhi.m
%>
%> @brief Assembles a matrix containing integrals over interior edges of 
%>        products of two basis functions.
%===============================================================================
%>
%> @brief Assembles the matrix @f$\mathsf{{S}}@f$ containing integrals
%>        over edges of products of two basis functions.
%>
%> The matrices @f$\mathsf{{S}} \in \mathbb{R}^{KN\times KN}@f$
%> consist of two kinds of contributions: diagonal blocks and off-diagonal 
%> blocks. Diagonal blocks are defined as 
%> @f[
%> [\mathsf{{S}}]_{(k-1)N+i,(k-1)N+j} =
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_D}
%>  \frac{1}{|E_{kn}|} \int_{E_{kn}} \varphi_{ki} \varphi_{kj} \mathrm{d}s \,,
%> @f]
%> and off-diagonal blocks are defined as
%> @f[
%> [\mathsf{{S}}]_{(k^--1)N+i,(k^+-1)N+j} =
%>  -\frac{1}{|E_{k^-n^-}|} \int_{E_{k^-n^-}} 
%>   \varphi_{k^-i} \varphi_{k^+j} \mathrm{d}s \,.
%> @f]
%> Entries in off-diagonal blocks are only non-zero for pairs of elements
%> @f$T_{k^-}, T_{k^+}@f$ with @f$\partial T_{k^-} \cap T_{k^+} \ne\emptyset@f$.
%> Note that the local edge index @f$n^-@f$ is given implicitly, since 
%> @f$\partial T_{k^-} \cap T_{k^+} \ne\emptyset@f$ consist of exactly one
%> edge @f$E_{k^-n^-} = E_{k^+n^+}@f$.
%> To allow for vectorization, the assembly is reformulated as
%> @f$\mathsf{{S}} = \mathsf{{S}}^\mathrm{diag} + 
%>    \mathsf{{S}}^\mathrm{offdiag}@f$ with the blocks defined as
%> @f[
%> \mathsf{{S}}^\mathrm{diag} = \sum_{n=1}^4
%>   \begin{bmatrix}
%>     \delta_{E_{1n}\in\mathcal{E}_\Omega} &   & \\
%>     & ~\ddots~ & \\
%>     &          & \delta_{E_{Kn}\in\mathcal{E}_\Omega}
%>   \end{bmatrix} \otimes [\hat{\mathsf{{S}}}^\mathrm{diag}]_{:,:,n}\;,
%> @f]
%> and
%> @f[
%> \mathsf{{S}}^\mathrm{offdiag} = -\sum_{n^-=1}^4
%>   \begin{bmatrix}
%>     0&\delta_{E_{1n^-} = E_{1n^+}}&\dots&\dots&\delta_{E_{1n^-}=E_{Kn^+}} \\
%>     \delta_{E_{2n^-} = E_{1n^+}}&0&\ddots& &\vdots \\
%>     \vdots & \ddots & \ddots & \ddots & \vdots \\
%>     \vdots & & \ddots & 0 & \delta_{E_{(K-1)n^-}=E_{Kn^+}} \\
%>     \delta_{E_{Kn^-} = E_{1n^+}}&\dots&\dots&\delta_{E_{Kn^-} = E_{(K-1)n^+}} &0
%>   \end{bmatrix} \otimes 
%> [\hat{\mathsf{{S}}}^\mathrm{offdiag}]_{:,:,n^-}\;,
%> @f]
%> where @f$\delta_{E_{kn}\in\mathcal{E}_\Omega}, \delta_{E_{in^-}=E_{jn^+}}@f$ 
%> denote the Kronecker delta and @f$\otimes@f$ denotes the Kronecker product.
%> The index @f$n^+@f$ is given implicitly by @f$n^-@f$ as described in
%> <code>mapLocalEdgeTetra()</code>.
%>
%> The entries of matrix 
%> @f$\hat{\mathsf{{S}}}^\mathrm{diag}\in\mathbb{R}^{N\times N\times4}@f$
%> are given by
%> @f[
%> [\hat{\mathsf{{S}}}^\mathrm{diag}]_{i,j,n} =
%>   \int_0^1 \hat{\varphi}_i \circ \hat{\mathbf{\gamma}}_n(s) 
%>   \hat{\varphi}_j\circ \hat{\mathbf{\gamma}}_n(s) \mathrm{d}s \,,
%> @f]
%> where the mapping @f$\hat{\mathbf{\gamma}}_n@f$ is defined in 
%> <code>gammaMapTetra()</code>. The entries of matrix
%> @f$\hat{\mathsf{{S}}}^\mathrm{offdiag} \in 
%>    \mathbb{R}^{N\times N\times 4}@f$ are defined as
%> @f[
%> [\hat{\mathsf{{S}}}^\mathrm{offdiag}]_{i,j,n^-} =
%>   \int_0^1 \hat{\varphi}_i \circ \hat{\mathbf{\gamma}}_{n^-}(s) 
%>   \hat{\varphi}_j\circ \hat{\mathbf{\gamma}}_{n^+}(s) \mathrm{d}s \,.
%> @f]
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>domainRectTrap()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0T    <code>logical</code> arrays that mark each elements
%>                    edges on which the matrix blocks should be
%>                    assembled @f$[K \times 4]@f$
%> @param refEdgePhiIntPhiInt  Local matrix 
%>                    @f$\hat{\mathsf{{S}}}^\text{diag}@f$ as provided
%>                    by <code>integrateRefEdgeTetraPhiIntPhiInt()</code>.
%>                    @f$[N \times N \times 4]@f$
%> @param refEdgePhiIntPhiExt Local matrix 
%>                    @f$\hat{\mathsf{{S}}}^\text{offdiag}@f$ as provided
%>                    by <code>integrateRefEdgeTetraPhiIntPhiExt()</code>.
%>                    @f$[N \times N \times 4]@f$
%> @param coefE0T     (optional) Coefficient vector that is applied to each
%>                    diagonal block. Defaults to <code>g.areaE0T</code>
%>                    @f$[K \times 4]@f$
%> @retval ret        The assembled matrix @f$[KN \times KN]@f$
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
function ret = assembleMatEdgeTetraPhiPhi(g, markE0T, refEdgePhiIntPhiInt, refEdgePhiIntPhiExt, coefE0T)
if nargin < 5
  coefE0T = g.areaE0T;
end % if
K = g.numT; N = size(refEdgePhiIntPhiInt, 1);
ret = sparse(K*N, K*N);
for n = 1 : 4
  areaE0Tint = markE0T(:,n) .* coefE0T(:,n);
  ret = ret + ...
        kron(spdiags(areaE0Tint, 0, K, K), refEdgePhiIntPhiInt(:,:,n)) - ...
        kron(bsxfun(@times, g.markE0TE0T{n}, areaE0Tint), refEdgePhiIntPhiExt(:,:,n));
end % for n
end % function