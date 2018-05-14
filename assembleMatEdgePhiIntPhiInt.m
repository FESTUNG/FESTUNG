% Assembles a matrix containing integrals over edges of products of two 
% basis functions from the interior of each element.

%===============================================================================
%> @file
%>
%> @brief Assembles a matrix containing integrals over edges of products of
%>        two basis functions from the interior of each element.
%===============================================================================
%>
%> @brief Assembles the matrix @f$\mathsf{{S}}_\mathrm{D}@f$ containing 
%>        integrals over edges of products of two basis functions from the 
%>        interior of each element.
%>
%> The matrix @f$\mathsf{{S}}_\mathrm{D} \in \mathbb{R}^{KN\times KN}@f$
%> is block diagonal and defined as 
%> @f[
%> [\mathsf{{S}}_\mathrm{D}]_{(k-1)N+i,(k-1)N+j} =
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_D}
%>  \int_{E_{kn}} \varphi_{ki} \varphi_{kj} \mathrm{d}s \,.
%> @f]
%> All other entries are zero.
%> To allow for vectorization, the assembly is reformulated as
%> @f[
%> \mathsf{{S}}_\mathrm{D} = \sum_{n=1}^{n_\mathrm{edges}}
%>   \begin{bmatrix}
%>     \delta_{E_{1n}\in\mathcal{E}_\mathrm{D}} &   & \\
%>     & ~\ddots~ & \\
%>     &          & \delta_{E_{Kn}\in\mathcal{E}_\mathrm{D}}
%>   \end{bmatrix} \otimes [\hat{\mathsf{{S}}}^\mathrm{diag}]_{:,:,n}\;,
%> @f]
%> where @f$\delta_{E_{kn}\in\mathcal{E}_\mathrm{D}}@f$ denotes the Kronecker 
%> delta and @f$\otimes@f$ denotes the Kronecker product.
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
%> <code>gammaMap()</code>.
%>
%> It is essentially the same as the diagonal part of
%> <code>assembleMatEdgePhiPhi()</code>.
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see, e.g., 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0T    <code>logical</code> arrays that mark each triangles
%>                    (boundary) edges on which the matrix blocks should be
%>                    assembled @f$[K \times n_\mathrm{edges}]@f$
%> @param refEdgePhiIntPhiInt  Local matrix 
%>                    @f$\hat{\mathsf{{S}}}^\text{diag}@f$ as provided
%>                    by <code>integrateRefEdgePhiIntPhiInt()</code> or
%>                    by <code>integrateRefEdgeTetraPhiIntPhiInt()</code>.
%>                    @f$[N \times N \times n_\mathrm{edges}]@f$
%> @param coefE0T     (optional) Coefficient vector that is applied to each
%>                    block. Defaults to <code>g.areaE0T</code>
%>                    @f$[K \times n_\mathrm{edges}]@f$
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
function ret = assembleMatEdgePhiIntPhiInt(g, markE0T, refEdgePhiIntPhiInt, coefE0T)
if nargin < 4, coefE0T = g.areaE0T; end

K = g.numT;  N = size(refEdgePhiIntPhiInt, 1); nEdges = size(g.E0T, 2);

% Check function arguments that are directly used
validateattributes(markE0T, {'logical'}, {'size', [K nEdges]}, mfilename, 'markE0T');
validateattributes(refEdgePhiIntPhiInt, {'numeric'}, {'size', [N N nEdges]}, mfilename, 'refEdgePhiIntPhiInt');

% Assemble matrix
ret = sparse(K*N, K*N);
for n = 1 : nEdges
  ret = ret + kron(spdiags(markE0T(:, n) .* coefE0T(:, n), 0, K, K), refEdgePhiIntPhiInt(:, :, n));
end % for
end % function
