% Assembles two matrices containing integrals over edges of products of a 
% linear Lagrange basis function and a DG basis function from the interior 
% of each element with a component of the edge normal.

%===============================================================================
%> @file assembleMatEdgePhiLagrPhiIntNu.m
%>
%> @brief Assembles two matrices containing integrals over edges of 
%>				products of a linear Lagrange basis function and a DG basis 
%>				function from the interior of each element with a component of 
%>				the edge normal.
%===============================================================================
%>
%> @brief Assembles two matrices containing integrals over edges of 
%>				products of a linear Lagrange basis function and a DG basis 
%>				function from the interior of each element with a component of 
%>				the edge normal.
%>
%> The matrix @f$\mathsf{{Q}_L}^m_\mathrm{N} \in \mathbb{R}^{L\times KN}@f$
%> is defined as 
%> @f[
%> [\mathsf{{Q}_L}^m_\mathrm{N}]_{i,(k-1)N+j} =
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_N}
%>  \int_{E_{kn}} \nu_{kn}^m \varphi_{i}^L \varphi_{kj} \mathrm{d}s \,.
%> @f]
%> with @f$\nu_{kn}@f$ the @f$m@f$-th component (@f$m\in\{1,2\}@f$) of the edge
%> normal. \varphi_i^L is the piecewise linear continuous function whose
%> value in the i-th global grid vertex is one and zero in all others.
%> All other entries are zero.
%> To allow for vectorization, the assembly is reformulated as
%> @f[
%> \mathsf{{Q}_L}^m_\mathrm{N} = \sum_{n=1}^3
%>   \begin{bmatrix}
%>     \delta_{E_{1n}\in\mathcal{E}_\mathrm{N}} &   & \\
%>     & ~\ddots~ & \\
%>     &          & \delta_{E_{Kn}\in\mathcal{E}_\mathrm{N}}
%>   \end{bmatrix} \circ \begin{bmatrix}
%>     \nu^m_{1n} | E_{1n} | &   & \\
%>     & ~\ddots~ & \\
%>     &          & \nu^m_{Kn} | E_{Kn} |
%>   \end{bmatrix} \otimes [\hat{\mathsf{{S}}}^\mathrm{diag}]_{:,:,n}\;,
%> @f]
%> where @f$\delta_{E_{kn}\in\mathcal{E}_\mathrm{N}}@f$ denotes the Kronecker 
%> delta, @f$\circ@f$ denotes the Hadamard product, and @f$\otimes@f$ denotes 
%> the Kronecker product.
%>
%> The entries of matrix 
%> @f$\hat{\mathsf{{S}}}^\mathrm{diag}\in\mathbb{R}^{N\times N\times3}@f$
%> are given by
%> @f[
%> [\hat{\mathsf{{S}}}^\mathrm{diag}]_{i,j,n} =
%>   \int_0^1 \hat{\varphi}_i \circ \hat{\mathbf{\gamma}}_n(s) 
%>   \hat{\varphi}_j\circ \hat{\mathbf{\gamma}}_n(s) \mathrm{d}s \,,
%> @f]
%> where the mapping @f$\hat{\mathbf{\gamma}}_n@f$ is defined in 
%> <code>gammaMap()</code>.
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0Tbdr <code>logical</code> arrays that mark each triangles
%>                    (boundary) edges on which the matrix blocks should be
%>                    assembled @f$[K \times 3]@f$
%> @param refEdgePhiIntPhiInt  Local matrix 
%>                    @f$\hat{\mathsf{{S}}}^\text{diag}@f$ as provided
%>                    by <code>integrateRefEdgePhiIntPhiInt()</code>.
%>                    @f$[N \times N \times 3]@f$
%> @param areaNuE0Tbdr (optional) argument to provide precomputed values
%>                    for the products of <code>markE0Tbdr</code>,
%>                    <code>g.areaE0T</code>, and <code>g.nuE0T</code>
%>                    @f$[3 \times 2 \text{ cell}]@f$
%> @retval ret        The assembled matrices @f$[2 \times 1 \text{ cell}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Hennes Hajduk, Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function ret = assembleVecEdgePhiLagrFuncQ0E0TNu(g, funcQ0E0T, markE0Tbdr, refEdgePhiLagrPerQuad, areaNuE0Tbdr)
% Extract dimensions
K = g.numT;

% Check function arguments that are directly used
validateattributes(markE0Tbdr, {'logical'}, {'size', [K 3]}, mfilename, 'markE0Tbdr');
validateattributes(refEdgePhiLagrPerQuad, {'numeric'}, {'size', [3 3 NaN]}, mfilename, 'refEdgePhiIntPhiInt');

if nargin > 4
  ret = assembleVecEdgePhiLagrFuncQ0E0TNu_withAreaNuE0Tbdr(g, funcQ0E0T, refEdgePhiLagrPerQuad, areaNuE0Tbdr);
elseif isfield(g, 'areaNuE0T')
  ret = assembleVecEdgePhiLagrFuncQ0E0TNu_noAreaNuE0Tbdr_withAreaNuE0T(g, funcQ0E0T, markE0Tbdr, refEdgePhiLagrPerQuad);
else
  ret = assembleVecEdgePhiLagrFuncQ0E0TNu_noAreaNuE0Tbdr_noAreaNuE0T(g, funcQ0E0T, markE0Tbdr, refEdgePhiLagrPerQuad);
end % if
end % function
%
%===============================================================================
%> @brief Helper function for the case that assembleMatEdgePhiLagrPhiIntNu()
%> was called with a precomputed field areaNuE0Tbdr.
%
function ret = assembleVecEdgePhiLagrFuncQ0E0TNu_withAreaNuE0Tbdr(g, funcQ0E0T, refEdgePhiLagrPerQuad, areaNuE0Tbdr)
% Extract dimensions
K = g.numT;  L = g.numV;

% Assemble matrices
ret = cell(2,1); 
ret{1} = sparse(L, 1); 
ret{2} = sparse(L, 1);
for n = 1	:	3
	markV0E0T = sparse(bsxfun(@eq, (1:L), g.V0T(:,mod(n,3)+1)))';
	ret{1} = ret{1} + kron(bsxfun(@times, markV0E0T, areaNuE0Tbdr{n,1}'), refEdgePhiLagrPerQuad(mod(n,3)+1,n,:));
  ret{2} = ret{2} + kron(bsxfun(@times, markV0E0T, areaNuE0Tbdr{n,2}'), refEdgePhiLagrPerQuad(mod(n,3)+1,n,:));
  
  markV0E0T = sparse(bsxfun(@eq, (1:L), g.V0T(:,mod(n+1,3)+1)))';
	ret{1} = ret{1} + kron(bsxfun(@times, markV0E0T, areaNuE0Tbdr{n,1}'), refEdgePhiLagrPerQuad(mod(n+1,3)+1,n,:));
  ret{2} = ret{2} + kron(bsxfun(@times, markV0E0T, areaNuE0Tbdr{n,2}'), refEdgePhiLagrPerQuad(mod(n+1,3)+1,n,:));
end % for
end % function
%
%===============================================================================
%> @brief Helper function for the case that assembleMatEdgePhiLagrPhiIntNu()
%> was called without a precomputed field areaNuE0Tbdr and parameter g provides
%> a precomputed field areaNuE0T.
%
function ret = assembleVecEdgePhiLagrFuncQ0E0TNu_noAreaNuE0Tbdr_withAreaNuE0T(g, funcQ0E0T, markE0Tbdr, refEdgePhiLagrPerQuad)
% Extract dimensions
K = g.numT;  L = g.numV;

% Assemble matrices
ret = cell(2,1); 
ret{1} = sparse(L, 1); 
ret{2} = sparse(L, 1);
for n = 1 : 3
  markV0E0T = sparse(bsxfun(@eq, (1:L), g.V0T(:,mod(n,3)+1)))'; 
  ret{1} = ret{1} + kron(bsxfun(@times, markV0E0T, markE0Tbdr(:,n)' .* g.areaNuE0T{n,1}'), refEdgePhiLagrPerQuad(mod(n,3)+1,n,:));
  ret{2} = ret{2} + kron(bsxfun(@times, markV0E0T, markE0Tbdr(:,n)' .* g.areaNuE0T{n,2}'), refEdgePhiLagrPerQuad(mod(n,3)+1,n,:));
  
  markV0E0T = sparse(bsxfun(@eq, (1:L), g.V0T(:,mod(n+1,3)+1)))'; 
  ret{1} = ret{1} + kron(bsxfun(@times, markV0E0T, markE0Tbdr(:,n)' .* g.areaNuE0T{n,1}'), refEdgePhiLagrPerQuad(mod(n+1,3)+1,n,:));
  ret{2} = ret{2} + kron(bsxfun(@times, markV0E0T, markE0Tbdr(:,n)' .* g.areaNuE0T{n,2}'), refEdgePhiLagrPerQuad(mod(n+1,3)+1,n,:));
end % for
end % function
%
%===============================================================================
%> @brief Helper function for the case that assembleMatEdgePhiLagrPhiIntNu()
%> was called without a precomputed field areaNuE0Tbdr and parameter g provides
%> no precomputed field areaNuE0T.
%
function ret = assembleVecEdgePhiLagrFuncQ0E0TNu_noAreaNuE0Tbdr_noAreaNuE0T(g, funcQ0E0T, markE0Tbdr, refEdgePhiLagrPerQuad)
% Extract dimensions
K = g.numT;  L = g.numV;

% Assemble matrices
ret = cell(2,1); 
ret{1} = sparse(L, 1); 
ret{2} = sparse(L, 1);
for n = 1 : 3
  Qkn1 = sparse( bsxfun(@times, bsxfun(@eq, (1:L), g.V0T(:,mod(n,3)+1)))', markE0Tbdr(:,n) .* g.areaE0T(:,n) .* g.nuE0T(:,n,1) );
  Qkn2 = sparse( bsxfun(@times, bsxfun(@eq, (1:L), g.V0T(:,mod(n,3)+1)))', markE0Tbdr(:,n) .* g.areaE0T(:,n) .* g.nuE0T(:,n,2) );
  ret{1} = ret{1} + kron(Qkn1, refEdgePhiLagrPerQuad(mod(n,3)+1,n,:));
  ret{2} = ret{2} + kron(Qkn2, refEdgePhiLagrPerQuad(mod(n,3)+1,n,:));
  
  Qkn1 = sparse( bsxfun(@times, bsxfun(@eq, (1:L), g.V0T(:,mod(n+1,3)+1)))', markE0Tbdr(:,n) .* g.areaE0T(:,n) .* g.nuE0T(:,n,1) );
  Qkn2 = sparse( bsxfun(@times, bsxfun(@eq, (1:L), g.V0T(:,mod(n+1,3)+1)))', markE0Tbdr(:,n) .* g.areaE0T(:,n) .* g.nuE0T(:,n,2) );
  ret{1} = ret{1} + kron(Qkn1, refEdgePhiLagrPerQuad(mod(n+1,3)+1,n,:));
  ret{2} = ret{2} + kron(Qkn2, refEdgePhiLagrPerQuad(mod(n+1,3)+1,n,:));
end % for
end % function
