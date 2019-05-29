% Assembles 24 matrices containing evaluations of one basis function in 
% quadrature points for each of the three local edges of each element as well as
% each combination of local edge indices of neighbouring elements multiplied
% with a component of the edge normal and the corresponding quadrature weight.

%===============================================================================
%> @file
%>
%> @brief Assembles 24 matrices containing evaluations of one basis function in 
%>        quadrature points for each of the three local edges of each element as
%>        well as each combination of local edge indices of neighbouring 
%>        elements multiplied with a component of the edge normal and the 
%>        corresponding quadrature weight.
%===============================================================================
%>
%> @brief Assembles matrices @f$\mathsf{{Q}}^{m,n}, m\in\{1,2\}, 
%>        n\in\{1,2,3\}@f$ containing evaluations of one basis function in 
%>        quadrature points for each of the three local edges of each element 
%>        multiplied with a component of the edge normal and the 
%>        corresponding quadrature weight
%>        as well as matrices @f$\mathsf{{Q}}^{m,n^-,n^+}, m\in\{1,2\}, 
%>        n^-,n^+\in\{1,2,3\}@f$ containing the same values but only 
%>        corresponding to triangles that share one edge of the appropriate 
%>        indices.
%>
%> The matrix @f$\mathsf{{Q}}^{m,n} \in \mathbb{R}^{KN\times KR}@f$ (R is the 
%> number of quadrature points and weights.) is block diagonal and defined as 
%> @f[
%> [\mathsf{{Q}}^{m,n}]_{(k-1)N+i,(k-1)R+r} = \frac{1}{2} \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_{\Omega}}
%>  \nu_{kn}^m \varphi_{ki}(q^r_{kn}) w^r_{kn} \,.
%> @f]
%> where @f$\nu_{kn}^m@f$ the @f$m@f$-th component (@f$m\in\{1,2\}@f$) of the edge
%> normal and @f$q^r_{kn}, w^r_{kn}@f$ the quadrature points and weights of edge @f$n@f$ of element @f$k@f$.
%> The matrix @f$\mathsf{{Q}}^{m,n^-,n^+} \in \mathbb{R}^{KN\times KR}@f$ is defined as
%> @f[
%> [\mathsf{{Q}}^{m,n^-,n^+}]_{(k^--1)N+i,(k^+-1)R+r} =
%>   \frac{1}{2} \nu^m_{k^-n^-} 
%>   \varphi_{k^-i}(q^r_{k^-n^-}) w^r_{k^-n^-} \,.
%> @f]
%> where @f$\nu_{k^-n^-}^m@f$ the @f$m@f$-th component (@f$m\in\{1,2\}@f$) of the edge
%> normal and @f$q^r_{k^-n^-}, w^r_{k^-n^-}@f$ the quadrature points and weights of edge @f$n^-@f$ of element @f$k^-@f$.
%> Entries in off-diagonal blocks are only non-zero for pairs of triangles
%> @f$T_{k^-}, T_{k^+}@f$ with @f$\partial T_{k^-} \cap T_{k^+} \ne\emptyset@f$.
%> Note that the local edge index @f$n^-@f$ is given implicitly, since 
%> @f$\partial T_{k^-} \cap T_{k^+} \ne\emptyset@f$ consist of exactly one
%> edge @f$E_{k^-n^-} = E_{k^+n^+}@f$.
%> All other entries are zero.
%> To allow for vectorization, the assembly is reformulated as
%> @f[
%> \mathsf{{Q}}^{m,n} = \frac{1}{2}
%>   \begin{bmatrix}
%>     \delta_{E_{1n}\in\mathcal{E}_{\Omega}} &   & \\
%>     & ~\ddots~ & \\
%>     &          & \delta_{E_{Kn}\in\mathcal{E}_{\Omega}}
%>   \end{bmatrix} \circ \begin{bmatrix}
%>     \nu^m_{1n} | E_{1n} | &   & \\
%>     & ~\ddots~ & \\
%>     &          & \nu^m_{Kn} | E_{Kn} |
%>   \end{bmatrix} 
%>  \otimes [\hat{\mathsf{{S}}}]_{:,:,n}\;,
%> @f]
%> and
%> @f[
%> \mathsf{{Q}}^{m,n^-,n^+} = \frac{1}{2} 
%>   \sum_{n^-=1}^3\sum_{n^+=1}^3
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
%> [\hat{\mathsf{{S}}}]_{:,:,n^-}\;,
%> @f]
%> where @f$\delta_{E_{kn}\in\mathcal{E}_{\Omega}}@f$ denotes the Kronecker 
%> delta, @f$\circ@f$ denotes the Hadamard product, and @f$\otimes@f$ denotes 
%> the Kronecker product.
%>
%> The entries of matrix 
%> @f$\hat{\mathsf{{S}}}\in\mathbb{R}^{N\times R \times 3}@f$
%> are given by
%> @f[
%> [\hat{\mathsf{{S}}}]_{i,r,n} =
%>   \hat{\varphi}_i \circ \hat{\mathbf{\gamma}}_n(\hat{q}^r) \hat{w}^r\,,
%> @f]
%> where the mapping @f$\hat{\mathbf{\gamma}}_n@f$ is defined in 
%> <code>gammaMap()</code> and @f$\hat{q}^r, \hat{w}^r@f$ are the 
%> quadrature points and weights of the interval @f$[0,1]@f$.
%>
%> @param g           The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param markE0Tint  <code>logical</code> arrays that mark each triangles
%>                    (boundary) edges on which the matrix blocks should be
%>                    assembled @f$[K \times 3]@f$
%> @param refEdgePhiIntPerQuad  Local matrix 
%>                    @f$\hat{\mathsf{S}}@f$ as provided
%>                    by <code>integrateRefEdgePhiIntPerQuad()</code>.
%>                    @f$[N \times R \times  3]@f$
%> @retval retDiag    The assembled matrices @f$[3 \times 2 \text{ cell}]@f$
%> @retval retOffdiag The assembled matrices @f$[3 \times 3 \times 2 \text{ cell}]@f$
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
function [retDiag, retOffdiag] = assembleMatEdgePhiNuPerQuad(g, markE0Tint, refEdgePhiIntPerQuad)

K = g.numT;

% Check function arguments that are directly used
validateattributes(markE0Tint, {'logical'}, {'size', [K 3]}, mfilename, 'markE0Tint');
validateattributes(refEdgePhiIntPerQuad, {'numeric'}, {'size', [NaN NaN 3]}, mfilename, 'refEdgePhiIntPerQuad');

retDiag = cell(3,2); retOffdiag = cell(3,3,2);
for nn = 1 : 3
	for np = 1 : 3
		if isfield(g, 'areaNuMarkE0TE0T')
			retOffdiag{nn,np,1} = 0.5*kron(g.areaNuMarkE0TE0T{nn,np,1}, refEdgePhiIntPerQuad(:,:,nn));
			retOffdiag{nn,np,2} = 0.5*kron(g.areaNuMarkE0TE0T{nn,np,2}, refEdgePhiIntPerQuad(:,:,nn));
		else
			retOffdiag{nn,np,1} = 0.5*kron(bsxfun(@times,g.markE0TE0T{nn,np},g.areaE0T(:,nn).*g.nuE0T(:,nn,1)), refEdgePhiIntPerQuad(:,:,nn));
			retOffdiag{nn,np,2} = 0.5*kron(bsxfun(@times,g.markE0TE0T{nn,np},g.areaE0T(:,nn).*g.nuE0T(:,nn,2)), refEdgePhiIntPerQuad(:,:,nn));
		end % if
	end % for
	if isfield(g, 'areaNuE0Tint')
		retDiag{nn,1} = 0.5*kron(spdiags(g.areaNuE0Tint{nn,1}, 0, K, K), refEdgePhiIntPerQuad(:,:,nn));
		retDiag{nn,2} = 0.5*kron(spdiags(g.areaNuE0Tint{nn,2}, 0, K, K), refEdgePhiIntPerQuad(:,:,nn));
	else
		retDiag{nn,1} = 0.5*kron(spdiags(g.areaE0T(:,nn).*g.nuE0T(:,nn,1).*markE0Tint(:,nn), 0, K, K), refEdgePhiIntPerQuad(:,:,nn));
		retDiag{nn,2} = 0.5*kron(spdiags(g.areaE0T(:,nn).*g.nuE0T(:,nn,2).*markE0Tint(:,nn), 0, K, K), refEdgePhiIntPerQuad(:,:,nn));
	end % if
end % for
end % function
