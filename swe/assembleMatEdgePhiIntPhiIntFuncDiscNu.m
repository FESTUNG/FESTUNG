% Assembles two matrices containing integrals over edges of products of two 
% basis functions and a discrete function from the interior of each element  
% with a component of the edge normal.

%===============================================================================
%> @file assembleMatEdgePhiIntPhiIntFuncDiscNu.m
%>
%> @brief Assembles two matrices containing integrals over edges of products of 
%>        two basis functions and a discrete function from the interior of each 
%>        element with a component of the edge normal.
%===============================================================================
%>
%> @brief Assembles two matrices @f$\mathsf{{Q}}^m_\mathrm{N}, m\in\{1,2\}@f$
%>        containing integrals over edges of products of two basis functions 
%>        and a discrete function from the interior of each element with a
%>        component of the edge normal.
%>
%> The matrix @f$\mathsf{{Q}}^m_\mathrm{N} \in \mathbb{R}^{KN\times KN}@f$
%> is block diagonal and defined as 
%> @f[
%> [\mathsf{{Q}}^m_\mathrm{N}]_{(k-1)N+i,(k-1)N+j} = sum_{l=1}^{N_\mathrm{data}}
%>  zbDisc(k,l) \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_N}
%>  \int_{E_{kn}} \nu_{kn}^m \varphi_{ki} \varphi_{kj} \varphi_{kl} \mathrm{d}s \,.
%> @f]
%> with @f$\nu_{kn}@f$ the @f$m@f$-th component (@f$m\in\{1,2\}@f$) of the edge
%> normal.
%> All other entries are zero.
%> To allow for vectorization, the assembly is reformulated as
%> @f[
%> \mathsf{{Q}}^m_\mathrm{N} = \sum_{n=1}^3 \sum_{l=1}^{N_\mathrm{data}}
%>   \begin{bmatrix}
%>     \delta_{E_{1n}\in\mathcal{E}_\mathrm{N}} &   & \\
%>     & ~\ddots~ & \\
%>     &          & \delta_{E_{Kn}\in\mathcal{E}_\mathrm{N}}
%>   \end{bmatrix} \circ \begin{bmatrix}
%>     \nu^m_{1n} | E_{1n} | &   & \\
%>     & ~\ddots~ & \\
%>     &          & \nu^m_{Kn} | E_{Kn} |
%>   \end{bmatrix} \circ
%>   \begin{bmatrix} zbDisc(1,l) &  & \\
%>   & ~\ddots~ & \\
%>   &          & zbDisc(K,l) \end{bmatrix}
%>  \otimes [\hat{\mathsf{{S}}}^\mathrm{diag}]_{:,:,l,n}\;,
%> @f]
%> where @f$\delta_{E_{kn}\in\mathcal{E}_\mathrm{N}}@f$ denotes the Kronecker 
%> delta, @f$\circ@f$ denotes the Hadamard product, and @f$\otimes@f$ denotes 
%> the Kronecker product.
%>
%> The entries of matrix 
%> @f$\hat{\mathsf{{S}}}^\mathrm{diag}\in\mathbb{R}^{N\times N \times {N_\mathrm{data}} \times 3}@f$
%> are given by
%> @f[
%> [\hat{\mathsf{{S}}}^\mathrm{diag}]_{i,j,l,n} =
%>   \int_0^1 \hat{\varphi}_i \circ \hat{\mathbf{\gamma}}_n(s) 
%>   \hat{\varphi}_j\circ \hat{\mathbf{\gamma}}_n(s) 
%>   \hat{\varphi}_l \circ \hat{\mathbf{\gamma}}_n(s) \mathrm{d}s \,,
%> @f]
%> where the mapping @f$\hat{\mathbf{\gamma}}_n@f$ is defined in 
%> <code>gammaMap()</code>.
%>
%> It is essentially the same as the diagonal part of
%> <code>assembleMatEdgePhiPhiPhiLinNu()</code>.
%>
%> @param g           The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param markE0Tbdr  <code>logical</code> arrays that mark each triangles
%>                    (boundary) edges on which the matrix blocks should be
%>                    assembled @f$[K \times 3]@f$
%> @param refEdgePhiIntPhiIntPhiInt  Local matrix 
%>                    @f$\hat{\mathsf{{S}}}^\text{diag}@f$ as provided
%>                    by <code>integrateRefEdgePhiIntPhiInt()</code>.
%>                    @f$[N \times N \times {N_\mathrm{data}}\times 3]@f$
%> @param dataDisc    A representation of the discrete function, e.g., as 
%>                    computed by <code>projectFuncCont2DataDisc()</code>
%> @param areaNuE0Tbdr (optional) argument to provide precomputed values
%>                    for the products of <code>markE0Tbdr</code>,
%>                    <code>g.areaE0T</code>, and <code>g.nuE0T</code>
%>                    @f$[3 \times 2 \text{ cell}]@f$
%> @retval ret        The assembled matrices @f$[2 \times 1 \text{ cell}]@f$
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
function ret = assembleMatEdgePhiIntPhiIntFuncDiscNu(g, markE0Tbdr, refEdgePhiIntPhiIntPhiInt, dataDisc, areaNuE0Tbdr)
% Extract dimensions
dataN = size(dataDisc, 2);
N = size(refEdgePhiIntPhiIntPhiInt, 1);

% Check function arguments that are directly used
validateattributes(markE0Tbdr, {'logical'}, {'size', [g.numT dataN]}, mfilename, 'markE0Tbdr');
validateattributes(refEdgePhiIntPhiIntPhiInt, {'numeric'}, {'size', [N N dataN 3]}, mfilename, 'refEdgePhiIntPhiIntPhiLin');

if nargin > 4
  ret = assembleMatEdgePhiIntPhiIntFuncNu_withAreaNuE0Tbdr(refEdgePhiIntPhiIntPhiInt, areaNuE0Tbdr, dataDisc);
elseif isfield(g, 'areaNuE0T')
  ret = assembleMatEdgePhiIntPhiIntFuncNu_noAreaNuE0Tbdr_withAreaNuE0T(g, markE0Tbdr, refEdgePhiIntPhiIntPhiInt, dataDisc);
else
  ret = assembleMatEdgePhiIntPhiIntFuncNu_noAreaNuE0Tbdr_noAreaNuE0T(g, markE0Tbdr, refEdgePhiIntPhiIntPhiInt, dataDisc);
end % if
end % function
%
%===============================================================================
%> @brief Helper function for the case that assembleMatEdgePhiIntPhiIntFuncDiscNu()
%> was called with a precomputed field areaNuE0Tbdr.
%
function ret = assembleMatEdgePhiIntPhiIntFuncNu_withAreaNuE0Tbdr(refEdgePhiIntPhiIntPhiInt, areaNuE0Tbdr, dataDisc)
% Extract dimensions
[K, dataN] = size(dataDisc);
N = size(refEdgePhiIntPhiIntPhiInt, 1);

% Assemble matrices
ret = cell(2,1); 
ret{1} = sparse(K*N, K*N); 
ret{2} = sparse(K*N, K*N);
for n = 1 : 3
  for l = 1 : dataN
    ret{1} = ret{1} + kron(spdiags(areaNuE0Tbdr{n,1} .* dataDisc(:,l), 0,K,K), refEdgePhiIntPhiIntPhiInt(:,:,l,n));
    ret{2} = ret{2} + kron(spdiags(areaNuE0Tbdr{n,2} .* dataDisc(:,l), 0,K,K), refEdgePhiIntPhiIntPhiInt(:,:,l,n));
  end % for
end % for
end % function
%
%===============================================================================
%> @brief Helper function for the case that assembleMatEdgePhiIntPhiIntFuncDiscNu()
%> was called without a precomputed field areaNuE0Tbdr and parameter g provides
%> a precomputed field areaNuE0T.
%
function ret = assembleMatEdgePhiIntPhiIntFuncNu_noAreaNuE0Tbdr_withAreaNuE0T(g, markE0Tbdr, refEdgePhiIntPhiIntPhiInt, dataDisc)
% Extract dimensions
[K, dataN] = size(dataDisc);
N = size(refEdgePhiIntPhiIntPhiInt, 1);

% Assemble matrices
ret = cell(2,1); 
ret{1} = sparse(K*N, K*N); 
ret{2} = sparse(K*N, K*N);
for n = 1 : 3
  for l = 1 : dataN
    ret{1} = ret{1} + kron(spdiags(markE0Tbdr(:,n).*g.areaNuE0T{n,1} .* dataDisc(:,l), 0,K,K), refEdgePhiIntPhiIntPhiInt(:,:,l,n));
    ret{2} = ret{2} + kron(spdiags(markE0Tbdr(:,n).*g.areaNuE0T{n,2} .* dataDisc(:,l), 0,K,K), refEdgePhiIntPhiIntPhiInt(:,:,l,n));
  end % for
end % for
end % function
%
%===============================================================================
%> @brief Helper function for the case that assembleMatEdgePhiIntPhiIntFuncDiscNu()
%> was called without a precomputed field areaNuE0Tbdr and parameter g provides
%> no precomputed field areaNuE0T.
%
function ret = assembleMatEdgePhiIntPhiIntFuncNu_noAreaNuE0Tbdr_noAreaNuE0T(g, markE0Tbdr, refEdgePhiIntPhiIntPhiInt, dataDisc)
% Extract dimensions
[K, dataN] = size(dataDisc);
N = size(refEdgePhiIntPhiIntPhiInt, 1);

% Assemble matrices
ret = cell(2,1); 
ret{1} = sparse(K*N, K*N); 
ret{2} = sparse(K*N, K*N);
for n = 1 : 3
  QNkn = markE0Tbdr(:,n) .* g.areaE0T(:,n);
  for l = 1 : dataN
    ret{1} = ret{1} + kron(spdiags(QNkn .* g.nuE0T(:,n,1) .* dataDisc(:,l), 0,K,K), refEdgePhiIntPhiIntPhiInt(:,:,l,n));
    ret{2} = ret{2} + kron(spdiags(QNkn .* g.nuE0T(:,n,2) .* dataDisc(:,l), 0,K,K), refEdgePhiIntPhiIntPhiInt(:,:,l,n));
  end % for
end % for
end % function
