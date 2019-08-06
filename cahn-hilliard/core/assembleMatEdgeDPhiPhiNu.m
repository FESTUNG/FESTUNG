% Assembles a matrix containing integrals over interior edges of products of a 
% derivative of a basis function, a basis function, and the normal of the triangle.

%===============================================================================
%> @file ./core/assembleMatEdgeDPhiPhiNu.m
%>
%> @brief Assembles a matrix containing integrals over interior edges of
%>        products of a derivative of a basis function, a basis function,
%>        and the normal of the triangle.
%===============================================================================
%>
%> @brief Assembles a matrix containing integrals over interior edges of
%>        products of a derivative of a basis function, a basis function,
%>        and the normal of the triangle.
%>
%> The matrices @f$\mathsf{{U}} \in \mathbb{R}^{KN\times KN}@f$
%> consist of two kinds of contributions: diagonal blocks and off-diagonal 
%> blocks. Diagonal blocks are defined as 
%> @f[
%> [\mathsf{{U}}]_{(k-1)N+i,(k-1)N+j} =
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_\Omega}
%>  \int_{E_{kn}} \nabla \varphi_{ki} \cdot \nu_{e_{kn}} \varphi_{kj} \mathrm{d}s \,,
%> @f]
%> and off-diagonal blocks are defined as
%> @f[
%> [\mathsf{{U}}]_{(k^--1)N+i,(k^+-1)N+j} =
%>  -\int_{E_{k^-n^-}} \nabla \varphi_{k^-i} \cdot \nu_{e_{k^-n^-}} \varphi_{k^+j} \mathrm{d}s \,.
%> @f]
%> Entries in off-diagonal blocks are potentially non-zero for pairs of elements
%> @f$T_{k^-}, T_{k^+}@f$ with @f$\partial T_{k^-} \cap T_{k^+} \ne\emptyset@f$.
%> Note that the local edge index @f$n^-@f$ is given implicitly, since 
%> @f$\partial T_{k^-} \cap T_{k^+} \ne\emptyset@f$ consist of exactly one
%> edge @f$E_{k^-n^-} = E_{k^+n^+}@f$; all possible combinations of @f$n^-@f$ and @f$n^+@f$
%> must be considered.
%>
%> To allow for vectorization, the assembly is reformulated as described in
%> [FrankKuzminRupp2018].
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see, e.g., 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0T    <code>logical</code> arrays that mark each elements
%>                    edges on which the matrix blocks should be
%>                    assembled @f$[K \times n_\mathrm{edges}]@f$
%> @param refEdgeDPhiIntPhiInt  Local matrix 
%>                    @f$\hat{\mathsf{{U}}}^\text{diag}@f$ as provided
%>                    by <code>integrateRefEdgeDPhiIntPhiInt()</code>.
%>                    @f$[N \times N \times n_\mathrm{edges}]@f$
%> @param refEdgeDPhiIntPhiExt Local matrix 
%>                    @f$\hat{\mathsf{{U}}}^\text{offdiag}@f$ as provided
%>                    by <code>integrateRefDEdgePhiIntPhiExt()</code>.
%>                    @f$[N \times N \times n_\mathrm{edges} \times n_\mathrm{edges}]@f$
%> @retval ret        The assembled matrix @f$[KN \times KN]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Florian Frank, Balthasar Reuter, Vadym Aizinger
%>
%> @author Andreas Rupp, 2018.
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
function ret = assembleMatEdgeDPhiPhiNu(g, markE0T, refEdgeDPhiIntPhiInt, refEdgeDPhiIntPhiExt, coefE0T)
K = g.numT;  N = size(refEdgeDPhiIntPhiInt{1}, 1); nEdges = size(g.E0T, 2);
if nargin < 5
  coefE0T = ones(K,nEdges);
end % if

% Check function arguments that are directly used
validateattributes(markE0T, {'logical'}, {'size', [K nEdges]}, mfilename, 'markE0T');
% validateattributes(refEdgeDPhiIntPhiInt, {'numeric'}, {'size', [N N nEdges]}, mfilename, 'refEdgePhiIntPhiInt');

ret = sparse(K*N, K*N);

if numel(g.markE0TE0T) == nEdges % mapping from nn to np implicitly given

  assert(0 == 1, 'Not yet implemented!');
  
else % mapping from nn to np explicitly given
  
%   validateattributes(refEdgeDPhiIntPhiExt, {'numeric'}, {'size', [N N nEdges nEdges 2]}, mfilename, 'refEdgePhiIntPhiExt');

  for nn = 1 : nEdges
    markCoefE0T = markE0T(:,nn) .* g.areaE0T(:,nn) .* coefE0T(:,nn) ./ g.detJ0T;
    ret = ret ...
            + kron(spdiags(markCoefE0T .* ( g.B(:,2,2) .* g.nuE0T(:,nn,1) - g.B(:,1,2) .* g.nuE0T(:,nn,2) ) , 0, K, K), refEdgeDPhiIntPhiInt{1}(:,:,nn)) ...
            + kron(spdiags(markCoefE0T .* (-g.B(:,2,1) .* g.nuE0T(:,nn,1) + g.B(:,1,1) .* g.nuE0T(:,nn,2) ) , 0, K, K), refEdgeDPhiIntPhiInt{2}(:,:,nn));
    for np = 1 : nEdges
      ret = ret ...
            - kron(bsxfun(@times, g.markE0TE0T{nn, np}, markCoefE0T .* ( g.B(:,2,2) .* g.nuE0T(:,nn,1) - g.B(:,1,2) .* g.nuE0T(:,nn,2) )), refEdgeDPhiIntPhiExt{1}(:,:,nn,np)) ...
            - kron(bsxfun(@times, g.markE0TE0T{nn, np}, markCoefE0T .* (-g.B(:,2,1) .* g.nuE0T(:,nn,1) + g.B(:,1,1) .* g.nuE0T(:,nn,2) )), refEdgeDPhiIntPhiExt{2}(:,:,nn,np));
    end % for np
  end % for nn
  
end % if
end % function
