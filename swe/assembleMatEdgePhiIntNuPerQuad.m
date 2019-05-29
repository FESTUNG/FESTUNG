% Assembles six matrices containing evaluations of one basis function in 
% quadrature points for each of the three local edges of each element multiplied
% with a component of the edge normal and the corresponding quadrature weight.

%===============================================================================
%> @file
%>
%> @brief Assembles six matrices containing evaluations of one basis function in 
%>        quadrature points for each of the three local edges of each element 
%>        multiplied with a component of the edge normal and the corresponding 
%>        quadrature weight.
%===============================================================================
%>
%> @brief Assembles matrices @f$\mathsf{{Q}}^{m,n}_\mathrm{L}, m\in\{1,2\}, 
%>        n\in\{1,2,3\}@f$ containing evaluations of one basis function in 
%>        quadrature points for each of the three local edges of each element 
%>        multiplied with a component of the edge normal and the corresponding 
%>        quadrature weight.
%>
%> The matrix @f$\mathsf{{Q}}^{m,n}_\mathrm{L} \in \mathbb{R}^{KN\times KR}@f$ (R is the 
%> number of quadrature points and weights.) is block diagonal and defined as 
%> @f[
%> [\mathsf{{Q}}^{m,n}_\mathrm{L}]_{(k-1)N+i,(k-1)R+r} = \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_L}
%>  \nu_{kn}^m \varphi_{ki}(q_{kn}^r) w_{kn}^r \,.
%> @f]
%> with @f$\nu_{kn}^m@f$ the @f$m@f$-th component (@f$m\in\{1,2\}@f$) of the edge
%> normal and @f$q^r_{kn},~w^r_{kn}@f$ the quadrature points and weights of edge @f$n@f$ of element @f$k@f$.
%> All other entries are zero.
%> To allow for vectorization, the assembly is reformulated as
%> @f[
%> \mathsf{{Q}}^{m,n}_\mathrm{L} = 
%>   \begin{bmatrix}
%>     \delta_{E_{1n}\in\mathcal{E}_\mathrm{L}} &   & \\
%>     & ~\ddots~ & \\
%>     &          & \delta_{E_{Kn}\in\mathcal{E}_\mathrm{L}}
%>   \end{bmatrix} \circ \begin{bmatrix}
%>     \nu^m_{1n} | E_{1n} | &   & \\
%>     & ~\ddots~ & \\
%>     &          & \nu^m_{Kn} | E_{Kn} |
%>   \end{bmatrix} 
%>  \otimes [\hat{\mathsf{{S}}}]_{:,:,n}\;,
%> @f]
%> where @f$\delta_{E_{kn}\in\mathcal{E}_\mathrm{L}}@f$ denotes the Kronecker 
%> delta, @f$\circ@f$ denotes the Hadamard product, and @f$\otimes@f$ denotes 
%> the Kronecker product.
%>
%> The entries of matrix 
%> @f$\hat{\mathsf{{S}}}\in\mathbb{R}^{N\times R \times 3}@f$
%> are given by
%> @f[
%> [\hat{\mathsf{{S}}}]_{i,r,n} =
%>   \hat{\varphi}_i \circ \hat{\mathbf{\gamma}}_n(\hat{q}^r) \hat{w}^r \,,
%> @f]
%> where the mapping @f$\hat{\mathbf{\gamma}}_n@f$ is defined in 
%> <code>gammaMap()</code> and \hat{q}^r, \hat{w}^r are the quadrature points 
%> and weights of edge of the interval @f$[0,1]@f$.
%>
%> It is essentially the same as the diagonal part of
%> <code>assembleMatEdgePhiNuPerQuad()</code>.
%>
%> @param g                     The lists describing the geometric and topological 
%>                              properties of a triangulation (see 
%>                              <code>generateGridData()</code>) 
%>                              @f$[1 \times 1 \text{ struct}]@f$
%> @param markE0Tbdr            <code>logical</code> arrays that mark each triangles
%>                              (boundary) edges on which the matrix blocks should be
%>                              assembled @f$[K \times 3]@f$
%> @param refEdgePhiIntPerQuad  Local matrix 
%>                              @f$\hat{\mathsf{S}}@f$ as provided
%>                              by <code>integrateRefEdgePhiIntPerQuad()</code>.
%>                              @f$[N \times R \times  3]@f$
%> @param areaNuE0Tbdr          (optional) argument to provide precomputed values
%>                              for the products of <code>markE0Tbdr</code>,
%>                              <code>g.areaE0T</code>, and <code>g.nuE0T</code>
%>                              @f$[3 \times 2 \text{ cell}]@f$
%> @retval ret                  The assembled matrices @f$[3 \times 2 \text{ cell}]@f$
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
function ret = assembleMatEdgePhiIntNuPerQuad(g, markE0Tbdr, refEdgePhiIntPerQuad, areaNuE0Tbdr)

K = g.numT;

% Check function arguments that are directly used
validateattributes(markE0Tbdr, {'logical'}, {'size', [K 3]}, mfilename, 'markE0Tbdr');
validateattributes(refEdgePhiIntPerQuad, {'numeric'}, {'size', [NaN NaN 3]}, mfilename, 'refEdgePhiIntPerQuad');

ret = cell(3,2);
for n = 1:3
	if nargin > 3
		ret{n,1} = kron(spdiags(areaNuE0Tbdr{n,1}, 0, K, K), refEdgePhiIntPerQuad(:,:,n));
		ret{n,2} = kron(spdiags(areaNuE0Tbdr{n,2}, 0, K, K), refEdgePhiIntPerQuad(:,:,n));
	else
		if isfield(g, 'areaNuE0T')
			ret{n,1} = kron(spdiags(g.areaNuE0T{n,1} .* markE0Tbdr(:,n), 0, K, K), refEdgePhiIntPerQuad(:,:,n));
			ret{n,2} = kron(spdiags(g.areaNuE0T{n,2} .* markE0Tbdr(:,n), 0, K, K), refEdgePhiIntPerQuad(:,:,n));
		else
			ret{n,1} = kron(spdiags(g.areaE0T(:,n) .* g.nuE0T(:,n,1) .* markE0Tbdr(:,n), 0, K, K), refEdgePhiIntPerQuad(:,:,n));
			ret{n,2} = kron(spdiags(g.areaE0T(:,n) .* g.nuE0T(:,n,2) .* markE0Tbdr(:,n), 0, K, K), refEdgePhiIntPerQuad(:,:,n));
		end % if
	end % if
end % for
end % function
