% Assembles a matrix containing integrals over interior edges of products of two 
% basis functions.

%===============================================================================
%> @file assembleMatEdgePhiPhi.m
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
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_\Omega}
%>  \int_{E_{kn}} \varphi_{ki} \varphi_{kj} \mathrm{d}s \,,
%> @f]
%> and off-diagonal blocks are defined as
%> @f[
%> [\mathsf{{S}}]_{(k^--1)N+i,(k^+-1)N+j} =
%>  -\int_{E_{k^-n^-}} \varphi_{k^-i} \varphi_{k^+j} \mathrm{d}s \,.
%> @f]
%> Entries in off-diagonal blocks are potentially non-zero for pairs of elements
%> @f$T_{k^-}, T_{k^+}@f$ with @f$\partial T_{k^-} \cap T_{k^+} \ne\emptyset@f$.
%> Note that the local edge index @f$n^-@f$ is given implicitly, since 
%> @f$\partial T_{k^-} \cap T_{k^+} \ne\emptyset@f$ consist of exactly one
%> edge @f$E_{k^-n^-} = E_{k^+n^+}@f$.
%> Depending on the element type, @f$n^+@f$ can either be deduced directly from
%> @f$n^-@f$ alone or all possible combinations of @f$n^-@f$ and @f$n^+@f$
%> must be considered.
%>
%> To allow for vectorization, the assembly is reformulated as
%> @f$\mathsf{{S}} = \mathsf{{S}}^\mathrm{diag} + \mathsf{{S}}^\mathrm{offdiag}@f$ 
%> with the blocks defined as
%> @f[
%> \mathsf{{S}}^\mathrm{diag} = \sum_{n=1}^{n_\mathrm{edges}}
%>   \begin{bmatrix}
%>     \delta_{E_{1n}\in\mathcal{E}_\Omega} &   & \\
%>     & ~\ddots~ & \\
%>     &          & \delta_{E_{Kn}\in\mathcal{E}_\Omega}
%>   \end{bmatrix} \otimes [\hat{\mathsf{{S}}}^\mathrm{diag}]_{:,:,n}\;,
%> @f]
%> and
%> @f[
%> \mathsf{{S}}^\mathrm{offdiag} = 
%>   -\sum_{n^-=1}^{n_\mathrm{edges}}\sum_{n^+=1}^{n_\mathrm{edges}}
%>   \begin{bmatrix}
%>     0&\delta_{E_{1n^-} = E_{2n^+}}&\dots&\dots&\delta_{E_{1n^-}=E_{Kn^+}} \\
%>     \delta_{E_{2n^-} = E_{1n^+}}&0&\ddots& &\vdots \\
%>     \vdots & \ddots & \ddots & \ddots & \vdots \\
%>     \vdots & & \ddots & 0 & \delta_{E_{(K-1)n^-}=E_{Kn^+}} \\
%>     \delta_{E_{Kn^-} = E_{1n^+}}&\dots&\dots&\delta_{E_{Kn^-} = E_{(K-1)n^+}} &0
%>   \end{bmatrix} \otimes 
%> [\hat{\mathsf{{S}}}^\mathrm{offdiag}]_{:,:,n^-,n^+}\;,
%> @f]
%> where @f$\delta_{E_{kn}\in\mathcal{E}_\Omega}, \delta_{E_{in^-}=E_{jn^+}}@f$ 
%> denote the Kronecker delta and @f$\otimes@f$ denotes the Kronecker product.
%>
%> The entries of matrix 
%> @f$\hat{\mathsf{{S}}}^\mathrm{diag}\in\mathbb{R}^{N\times N\times{n_\mathrm{edges}}}@f$
%> are given by
%> @f[
%> [\hat{\mathsf{{S}}}^\mathrm{diag}]_{i,j,n} =
%>   \int_0^1 \hat{\varphi}_i \circ \hat{\mathbf{\gamma}}_n(s) 
%>   \hat{\varphi}_j\circ \hat{\mathbf{\gamma}}_n(s) \mathrm{d}s \,,
%> @f]
%> where the mapping @f$\hat{\mathbf{\gamma}}_n@f$ is defined in 
%> <code>gammaMap()</code>. The entries of matrix
%> @f$\hat{\mathsf{{S}}}^\mathrm{offdiag} \in 
%>    \mathbb{R}^{N\times N\times {n_\mathrm{edges}}\times {n_\mathrm{edges}}}@f$ 
%> are defined as
%> @f[
%> [\hat{\mathsf{{S}}}^\mathrm{offdiag}]_{i,j,n^-,n^+} =
%>   \int_0^1 \hat{\varphi}_i \circ \hat{\mathbf{\gamma}}_{n^-}(s) 
%>   \hat{\varphi}_j\circ \hat{\mathbf{\vartheta}}_{n^-n^+} \circ
%>   \hat{\mathbf{\gamma}}_{n^-}(s) \mathrm{d}s \,,
%> @f]
%> with the mapping @f$\hat{\mathbf{\vartheta}}_{n^-n^+}@f$ as described in
%> <code>theta()</code>.
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see, e.g., 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0T    <code>logical</code> arrays that mark each elements
%>                    edges on which the matrix blocks should be
%>                    assembled @f$[K \times n_\mathrm{edges}]@f$
%> @param refEdgePhiIntPhiInt  Local matrix 
%>                    @f$\hat{\mathsf{{S}}}^\text{diag}@f$ as provided
%>                    by <code>integrateRefEdgePhiIntPhiInt()</code>.
%>                    @f$[N \times N \times n_\mathrm{edges}]@f$
%> @param refEdgePhiIntPhiExt Local matrix 
%>                    @f$\hat{\mathsf{{S}}}^\text{offdiag}@f$ as provided
%>                    by <code>integrateRefEdgePhiIntPhiExt()</code>.
%>                    @f$[N \times N \times n_\mathrm{edges}]@f$ or
%>                    @f$[N \times N \times n_\mathrm{edges} \times n_\mathrm{edges}]@f$
%> @retval ret        The assembled matrix @f$[KN \times KN]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Florian Frank, Balthasar Reuter, Vadym Aizinger
%>
%> @author Balthasar Reuter, 2017.
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
function ret = assembleMatEdgePhiPhi(g, markE0T, refEdgePhiIntPhiInt, refEdgePhiIntPhiExt, coefE0T)
if nargin < 5
  coefE0T = g.areaE0T;
end % if
K = g.numT;  N = size(refEdgePhiIntPhiInt, 1); nEdges = size(g.E0T, 2);

% Check function arguments that are directly used
validateattributes(markE0T, {'logical'}, {'size', [K nEdges]}, mfilename, 'markE0T');
validateattributes(refEdgePhiIntPhiInt, {'numeric'}, {'size', [N N nEdges]}, mfilename, 'refEdgePhiIntPhiInt');

ret = sparse(K*N, K*N);

if numel(g.markE0TE0T) == nEdges % mapping from nn to np implicitly given
  
  validateattributes(refEdgePhiIntPhiExt, {'numeric'}, {'size', [N N nEdges]}, mfilename, 'refEdgePhiIntPhiExt');
  
  for n = 1 : nEdges
    markCoefE0T = markE0T(:,n) .* coefE0T(:,n);
    ret = ret + ...
          kron(spdiags(markCoefE0T, 0, K, K), refEdgePhiIntPhiInt(:,:,n)) - ...
          kron(bsxfun(@times, g.markE0TE0T{n}, markCoefE0T), refEdgePhiIntPhiExt(:,:,n));
  end % for n
  
else % mapping from nn to np explicitly given
  
  validateattributes(refEdgePhiIntPhiExt, {'numeric'}, {'size', [N N nEdges nEdges]}, mfilename, 'refEdgePhiIntPhiExt');

  for nn = 1 : nEdges
    markCoefE0T = markE0T(:,nn) .* coefE0T(:,nn);
    ret = ret + kron(spdiags(markCoefE0T, 0, K, K), refEdgePhiIntPhiInt(:,:,nn));
    for np = 1 : nEdges
      ret = ret - kron(bsxfun(@times, g.markE0TE0T{nn, np}, markCoefE0T), refEdgePhiIntPhiExt(:,:,nn,np));
    end % for np
  end % for nn
  
end % if
end % function
