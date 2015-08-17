% Assembles a matrix containing integrals over edges of products of two 
% basis functions from the interior of each element.
%
%> @file assembleMatEdgePhiIntPhiInt.m
%>
%> @brief Assembles a matrix containing integrals over edges of products of
%>        two basis functions from the interior of each element.
%>
%> <code>ret = assembleMatEdgePhiIntPhiInt(g, markE0Tbdr, refEdgePhiIntPhiInt)
%> </code> assembles the matrix @f$\mathbf{\mathsf{S}}_\mathrm{D}@f$
%> containing integrals over edges of products of two basis functions from the 
%> interior of each element.
%>
%> The matrix @f$\mathbf{\mathsf{S}}_\mathrm{D}@f$ is defined as follows
%> @f[
%> \mathbf{\mathsf{S}}_\mathrm{D} = \sum_{n=1}^3
%>   \begin{bmatrix}
%>     \delta_{E_{1n}\in\mathcal{E}_\mathrm{D}} &   & \\
%>     & ~\ddots~ & \\
%>     &          & \delta_{E_{Kn}\in\mathcal{E}_\mathrm{D}}
%>   \end{bmatrix} \otimes [\hat{\mathbf{\mathsf{S}}}^\mathrm{diag}]_{:,:,n}\;,
%> @f]
%> where @f$\delta_{E_{kn}\in\mathcal{E}_\mathrm{D}}@f$ denotes the Kronecker 
%> delta and @f$\otimes@f$ denotes the Kronecker product.
%>
%> It is essentially the same as the diagonal part of
%> <code>assembleMatEdgePhiPhi()</code>.
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0Tbdr <code>logical</code> arrays that mark each triangles
%>                    (boundary) edges on which the matrix blocks should be
%>                    assembled @f$[K \times 3]@f$
%> @param refEdgePhiIntPhiInt  Local matrix 
%>                    @f$\hat{\mathbf{\mathsf{S}}}^\text{diag}@f$ as provided
%>                    by <code>integrateRefEdgePhiIntPhiInt()</code>
%>                    @f$[N \times N]@f$
%> @retval ret        The assembled matrix @f$[KN \times KN]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2015 Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function ret = assembleMatEdgePhiIntPhiInt(g, markE0Tbdr, refEdgePhiIntPhiInt)
K = g.numT;  N = size(refEdgePhiIntPhiInt, 1);
ret = sparse(K*N, K*N);
for n = 1 : 3
  ret = ret + kron(spdiags(markE0Tbdr(:,n),0,K,K), refEdgePhiIntPhiInt(:,:,n));
end % for
end % function
