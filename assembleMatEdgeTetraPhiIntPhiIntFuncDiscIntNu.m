% Assembles two matrices containing integrals over edges of products of two 
% basis functions from the interior of each element and a function in discrete
% representation with a component of the edge normal.

%===============================================================================
%> @file assembleMatEdgeTetraPhiIntPhiIntFuncDiscIntNu.m
%>
%> @brief Assembles two matrices containing integrals over edges of products of
%>        two basis functions from the interior of each element and a function
%>        in discrete representation with a component of the edge normal.
%===============================================================================
%>
%> @brief Assembles matrices
%>        @f$\mathsf{{R}}^m_\mathrm{D}, m\in\{1,2\}@f$ containing integrals over  
%>        edges of products of two basis functions and a function in discrete 
%>        representation from the interior of each element.
%>
%> The matrix @f$\mathsf{{R}}^m_\mathrm{D}\in\mathbb{R}^{KN\times KN}@f$
%> is block diagonal and defined as
%> @f[
%> [\mathsf{{R}}^m_\mathrm{D}]_{(k-1)N+i,(k-1)N+j} =
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_D}
%>  \nu_{kn}^m \sum_{l=1}^N D_{kl}(t) 
%>  \int_{E_kn} \varphi_{ki}\varphi_{kl}\varphi_{kj} \,,
%> @f]
%> with @f$\nu_{kn}@f$ the @f$m@f$-th component (@f$m\in\{1,2\}@f$) of the unit
%> normal and @f$D_{kl}(t)@f$ the coefficients of the discrete representation
%> of a function 
%> @f$ d_h(t, \mathbf{x}) = \sum_{l=1}^N D_{kl}(t) \varphi_{kl}(\mathbf{x}) @f$
%> on an element @f$T_k@f$.
%> All other entries are zero.
%> To allow for vectorization, the assembly is reformulated as
%> @f[
%> \mathsf{{R}}^m_\mathrm{D} = \sum_{n=1}^4 \sum_{l=1}^N
%>   \begin{bmatrix}
%>     \delta_{E_{1n}\in\mathcal{E}_\mathrm{D}} &   & \\
%>     & ~\ddots~ & \\
%>     &          & \delta_{E_{Kn}\in\mathcal{E}_\mathrm{D}}
%>   \end{bmatrix} \circ
%>   \begin{bmatrix}
%>     \nu^m_{1n} |E_{1n}| D_{1l}(t) & & \\
%>     & ~\ddots~ & \\
%>     &          & \nu^m_{Kn} |E_{Kn}| D_{Kl}(t)
%>   \end{bmatrix} \otimes [\hat{\mathsf{{R}}}^\mathrm{diag}]_{:,:,l,n}\;,
%> @f]
%> where @f$\delta_{E_{kn}\in\mathcal{E}_\mathrm{D}}@f$ denotes the Kronecker 
%> delta, @f$\circ@f$ denotes the Hadamard product, and @f$\otimes@f$ denotes
%> the Kronecker product.
%>
%> The entries of matrix 
%> @f$\hat{\mathsf{{R}}}^\mathrm{diag}\in\mathbb{R}^{N\times N\times N\times4}@f$
%> are given by
%> @f[
%> [\hat{\mathsf{{R}}}^\mathrm{diag}]_{i,j,l,n} =
%>   \int_0^1 \hat{\varphi}_i \circ \hat{\mathbf{\gamma}}_n(s) 
%>   \hat{\varphi}_l \circ \hat{\mathbf{\gamma}}_n(s) 
%>   \hat{\varphi}_j\circ \hat{\mathbf{\gamma}}_n(s) \mathrm{d}s \,,
%> @f]
%> where the mapping @f$\hat{\mathbf{\gamma}}_n@f$ is defined in 
%> <code>gammaMapTetra()</code>.
%>
%> It is essentially the same as the diagonal part of
%> <code>assembleMatEdgeTetraPhiPhiFuncDiscNu()</code>.
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>domainRectTrap()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0T    <code>logical</code> arrays that mark each elements
%>                    edges on which the matrix blocks should be
%>                    assembled @f$[K \times 4]@f$
%> @param refEdgePhiIntPhiIntPhiInt  Local matrix 
%>                    @f$\hat{\mathsf{{R}}}^\text{diag}@f$ as provided
%>                    by <code>integrateRefEdgeTetraPhiIntPhiIntPhiInt()</code>.
%>                    @f$[N \times N \times N \times 4]@f$
%> @param dataDisc    A representation of the discrete function 
%>                    @f$d_h(\mathbf(x))@f$, e.g., as computed by 
%>                    <code>projectFuncCont2DataDiscTetra()</code>
%>                    @f$[K \times N]@f$
%> @retval ret        The assembled matrices @f$[2 \times 1 \text{ cell}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
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
function ret = assembleMatEdgeTetraPhiIntPhiIntFuncDiscIntNu(g, markE0T, refEdgePhiIntPhiIntPhiInt, dataDisc)
if iscell(dataDisc)
  if isequal(size(dataDisc), [2 2])
    ret = assembleMatEdgeTrapPhiIntPhiIntFuncDiscMatrixIntNu(g, markE0T, refEdgePhiIntPhiIntPhiInt, dataDisc);
  elseif isvector(dataDisc) && numel(dataDisc) == 2
    ret = assembleMatEdgeTrapPhiIntPhiIntFuncDiscVectorIntNu(g, markE0T, refEdgePhiIntPhiIntPhiInt, dataDisc);
  else
    error('dataDisc must be a KxN-matrix or a 2x1-cell or a 2x2-cell of such matrices.')
  end % if
else
  ret = assembleMatEdgeTrapPhiIntPhiIntFuncDiscScalarIntNu(g, markE0T, refEdgePhiIntPhiIntPhiInt, dataDisc);
end % if
end % function

function ret = assembleMatEdgeTrapPhiIntPhiIntFuncDiscMatrixIntNu(g, markE0T, refEdgePhiIntPhiIntPhiInt, dataDisc)
[K, N] = size(dataDisc{1,1});
ret = { sparse(K*N, K*N), sparse(K*N, K*N) };
for n = 1 : 4
  for r = 1 : 2
    areaNuE0Tbdr = g.areaE0T(:,n) .* g.nuE0T(:,n,r) .* markE0T(:,n);
    for m = 1 : 2
      for l = 1 : N
        ret{m} = ret{m} + kron(spdiags(areaNuE0Tbdr .* dataDisc{r,m}(:,l), 0, K, K), ...
                               refEdgePhiIntPhiIntPhiInt(:,:,l,n));
      end % for l
    end  % for m
  end % for r
end  % for n
end  % function

function ret = assembleMatEdgeTrapPhiIntPhiIntFuncDiscVectorIntNu(g, markE0T, refEdgePhiIntPhiIntPhiInt, dataDisc)
[K, N] = size(dataDisc{1});
ret = { sparse(K*N, K*N), sparse(K*N, K*N) };
for n = 1 : 4
  for m = 1 : 2
    areaNuE0Tbdr = g.areaE0T(:,n) .* g.nuE0T(:,n,m) .* markE0T(:,n);
    for l = 1 : N
      ret{m} = ret{m} + kron(spdiags(areaNuE0Tbdr .* dataDisc{m}(:,l), 0, K, K), ...
                             refEdgePhiIntPhiIntPhiInt(:,:,l,n));
    end % for l
  end  % for m
end  % for n
end  % function

function ret = assembleMatEdgeTrapPhiIntPhiIntFuncDiscScalarIntNu(g, markE0T, refEdgePhiIntPhiIntPhiInt, dataDisc)
[K, N] = size(dataDisc);
ret = { sparse(K*N, K*N), sparse(K*N, K*N) };
for n = 1 : 4
  for m = 1 : 2
    areaNuE0Tbdr = g.areaE0T(:,n) .* g.nuE0T(:,n,m) .* markE0T(:,n);
    for l = 1 : N
      ret{m} = ret{m} + kron(spdiags(areaNuE0Tbdr .* dataDisc(:,l), 0, K, K), ...
                             refEdgePhiIntPhiIntPhiInt(:,:,l,n));
    end % for l
  end  % for m
end  % for n
end  % function