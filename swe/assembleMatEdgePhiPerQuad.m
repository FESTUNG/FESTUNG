% Assembles nine matrices containing evaluations of one basis function 
% in quadrature points multiplied with the corresponding quadrature 
% weight for each of the three local edges of each element which has a 
% neighbouring element to the local edge of the appropriate edge index.

%===============================================================================
%> @file
%>
%> @brief Assembles nine matrices containing evaluations of one basis function 
%>        in quadrature points multiplied with the corresponding quadrature 
%>        weight for each of the three local edges of each element which has a 
%>        neighbouring element to the local edge of the appropriate edge index.
%===============================================================================
%>
%> @brief Assembles matrices @f$\mathsf{{Q}}^{n^-,n^+},
%>        n^-,n^+\in\{1,2,3\}@f$ containing evaluations of one basis function 
%>        in quadrature points multiplied with the corresponding quadrature 
%>        weight for each of the three local edges of each element which has a 
%>        neighbouring element to the local edge of the appropriate edge index.
%>
%> The matrix @f$\mathsf{{Q}}^{n^-,n^+} \in \mathbb{R}^{KN\times KR}@f$ (R is the 
%> number of quadrature points and weights.) is block diagonal and defined as 
%> @f[
%> [\mathsf{{Q}}^{n^-n^+}]_{(k^--1)N+i,(k^--1)R+r} = \frac{1}{2} \left( \sum_{k^+=1}^K \delta_{E_{k^-n^-} = E_{k^+n^+}} \right)
%> \varphi_{k^-i}(q^r_{k^-n^-}) w^r_{k^-n^-} \,.
%> @f]
%> with @f$q^r_{k^-n^-}, w^r_{k^-n^-}@f$ the quadrature points and weights of edge @f$n^-@f$ of element @f$k^-@f$.
%>
%> All other entries are zero.
%> To allow for vectorization, the assembly is reformulated as
%> @f[
%> \mathsf{{Q}}^{n^-n^+} = \frac{1}{2}
%>   \begin{bmatrix}
%>     \sum_{k^+=1}^K \delta_{E_{1n^-} = E_{k^+n^+}} \\
%>     \vdots \\
%>     \sum_{k^+=1}^K \delta_{E_{Kn^-} = E_{k^+n^+}}
%>   \end{bmatrix} \circ \begin{bmatrix}
%>     | E_{1n} | &   & \\
%>     & ~\ddots~ & \\
%>     &          & | E_{Kn} |
%>   \end{bmatrix} 
%>  \otimes [\hat{\mathsf{{S}}}]_{:,:,n}\;,
%> @f]
%> where @f$\delta_{E_{k^-n^-} = E_{k^+n^+}}@f$ denotes the Kronecker 
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
%> @param refEdgePhiIntPerQuad  Local matrix 
%>                    @f$\hat{\mathsf{S}}@f$ as provided
%>                    by <code>integrateRefEdgePhiIntPerQuad()</code>.
%>                    @f$[N \times R \times  3]@f$
%> @retval retOffdiag The assembled matrices @f$[3 \times 3 \text{ cell}]@f$
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
%>
function ret = assembleMatEdgePhiPerQuad(g, refEdgePhiIntPerQuad)

% Check function arguments that are directly used
validateattributes(refEdgePhiIntPerQuad, {'numeric'}, {'size', [NaN NaN 3]}, mfilename, 'refEdgePhiIntPerQuad');

K = g.numT;
ret = cell(3,3);
for nn = 1 : 3
	for np = 1: 3
		if isfield(g, 'areaMarkE0T')
			ret{nn,np} = 0.5 * kron( spdiags(g.areaMarkE0T{nn,np}, 0, K, K), refEdgePhiIntPerQuad(:,:,nn) );
		else
			ret{nn,np} = 0.5 * kron( spdiags((g.markE0TE0T{nn,np} * ones(K,1)) .* g.areaE0T(:,nn), 0, K, K), refEdgePhiIntPerQuad(:,:,nn) );
		end % if
	end % for
end % for
end % function
