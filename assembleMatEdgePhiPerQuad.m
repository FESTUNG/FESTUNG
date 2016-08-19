% Assembles nine matrices containing evaluations of one basis function in 
% quadrature points for each of the three local edges of each element as well as
% each combination of local edge indices of neighbouring elements multiplied
% with the corresponding quadrature weight.

%===============================================================================
%> @file assembleMatEdgePhiPerQuad.m
%>
%> @brief Assembles nine matrices containing evaluations of one basis function 
%>        in quadrature points for each of the three local edges of each element
%>        as well as each combination of local edge indices of neighbouring 
%>        elements multiplied with the corresponding quadrature weight.
%===============================================================================
%>
%> @brief Assembles nine matrices containing evaluations of one basis function 
%>        in quadrature points for each of the three local edges of each element
%>        as well as each combination of local edge indices of neighbouring 
%>        elements multiplied with the corresponding quadrature weight.
%>
%> The matrix @f$\mathsf{{Q}}^m_{n^-n^+} \in \mathbb{R}^{KN\times KR}@f$ (R is the 
%> number of quadrature points and weights.) is block diagonal and defined as 
%> @f[
%> [\mathsf{{Q}}^m_{n^-n^+}]_{(k-1)N+i,(k-1)R+r} = \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_{\Omega}}
%> \varphi_{ki}(q^r_{kn}) w^r_{kn} \,.
%> @f]
%> with q^r_{kn}, w^r_{kn} the quadrature points and weights of edge n of element k.
%>
%> All other entries are zero.
%> To allow for vectorization, the assembly is reformulated as
%> @f[
%> \mathsf{{Q}}_{n^-n^+} =
%>   \\begin{bmatrix}
%>     sum_{k=1}^K 0&\delta_{E_{1n^-} = E_{kn^+}} \\
%>     \vdots \\
%>     sum_{k=1}^K 0&\delta_{E_{Kn^-} = E_{kn^+}}
%>   \end{bmatrix} \circ \begin{bmatrix}
%>     | E_{1n} | &   & \\
%>     & ~\ddots~ & \\
%>     &          & | E_{Kn} |
%>   \end{bmatrix} 
%>  \otimes [\hat{\mathsf{{S}}}]_{:,:,n}\;,
%> @f]
%> where @f$\delta_{E_{kn}\in\mathcal{E}_\mathrm{N}}@f$ denotes the Kronecker 
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
%> <code>gammaMap()</code> and \hat{q}^r, \hat{w}^r are the 
%> quadrature points and weights of edge @f$n@f$ of the reference element.
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
%
function ret = assembleMatEdgePhiPerQuad(g, refEdgePhiIntPerQuad)
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
