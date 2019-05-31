% Assembles two matrices containing integrals over interior edges of products of
% two basis functions with a discontinuous coefficient function and with a
% component of the edge normal.

%===============================================================================
%> @file
%>
%> @brief Assembles two matrices containing integrals over interior edges of 
%>        products of two basis functions with a discontinuous coefficient 
%>        function and with a component of the edge normal.
%===============================================================================
%>
%> @brief Assembles the matrices @f$\mathsf{{R}}^m, m\in\{1,2\}@f$ containing 
%>        integrals over interior edges of products of two basis functions with 
%>        a discontinuous coefficient function and a component of the edge normal.
%>
%> The matrices @f$\mathsf{{R}}^m\in\mathbb{R}^{KN\times KN}@f$
%> consist of two kinds of contributions: diagonal blocks and off-diagonal
%> blocks. Diagonal blocks are defined as
%> @f[
%> [\mathsf{{R}}^m]_{(k-1)N+i,(k-1)N+j} = \frac{1}{2} 
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_\Omega}
%>  \nu_{kn}^m \sum_{l=1}^N D_{kl}(t) 
%>  \int_{E_kn} \varphi_{ki}\varphi_{kl}\varphi_{kj} \,,
%> @f]
%> and off-diagonal blocks are defined as
%> @f[
%> [\mathsf{{R}}^m]_{(k^--1)N+i,(k^+-1)N+j} =
%>   \frac{1}{2} \nu^m_{k^-n^-} \sum_{l=1}^N D_{k^-l}(t) \int_{E_{k^-n^-}} 
%>   \varphi_{k^-i} \varphi_{k^-l} \varphi_{k^+j} \mathrm{d}s \,.
%> @f]
%> with @f$\nu_{kn}@f$ the @f$m@f$-th component (@f$m\in\{1,2\}@f$) of the unit
%> normal and @f$D_{kl}(t)@f$ the coefficients of the discrete representation
%> of a function 
%> @f$ d_h(t, \mathbf{x}) = \sum_{l=1}^N D_{kl}(t) \varphi_{kl}(\mathbf{x}) @f$
%> on a triangle @f$T_k@f$.
%> All other entries are zero.
%>
%> To allow for vectorization, the assembly is reformulated as
%> @f$\mathsf{{R}}^m = \mathsf{{R}}^{m,\mathrm{diag}} + 
%>    \mathsf{{R}}^{m,\mathrm{offdiag}}@f$ with the blocks defined as
%> @f[
%> \mathsf{{R}}^{m,\mathrm{diag}} = \frac{1}{2} \sum_{n=1}^3 \sum_{l=1}^N
%>   \begin{bmatrix}
%>     \delta_{E_{1n}\in\mathcal{E}_\Omega} &   & \\
%>     & ~\ddots~ & \\
%>     &          & \delta_{E_{Kn}\in\mathcal{E}_\Omega}
%>   \end{bmatrix} \circ
%>   \begin{bmatrix}
%>     \nu^m_{1n} |E_{1n}| D_{1l}(t) & & \\
%>     & ~\ddots~ & \\
%>     &          & \nu^m_{Kn} |E_{Kn}| D_{Kl}(t)
%>   \end{bmatrix} \otimes [\hat{\mathsf{{R}}}^\mathrm{diag}]_{:,:,l,n}\;,
%> @f]
%> and
%> @f[
%> \mathsf{{R}}^{m,\mathrm{offdiag}} = \frac{1}{2} 
%>   \sum_{n^-=1}^3\sum_{n^+=1}^3 \sum_{l=1}^N
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
%>   \end{bmatrix} \circ \begin{bmatrix}
%>     D_{1l}(t) & \dots & D_{1l}(t) \\
%>     \vdots & & \vdots \\
%>     D_{Kl}(t) & \dots & D_{Kl}(t)
%>   \end{bmatrix} \otimes 
%> [\hat{\mathsf{{R}}}^\mathrm{offdiag}]_{:,:,l,n^-,n^+}\;,
%> @f]
%> where @f$\delta_{E_{kn}\in\mathcal{E}_\Omega}, \delta_{E_{in^-}=E_{jn^+}}@f$ 
%> denote the Kronecker delta, @f$\circ@f$ denotes the Hadamard product, and 
%> @f$\otimes@f$ denotes the Kronecker product.
%>
%> The entries of matrix 
%> @f$\hat{\mathsf{{R}}}^\mathrm{diag}\in\mathbb{R}^{N\times N\times N\times3}@f$
%> are given by
%> @f[
%> [\hat{\mathsf{{R}}}^\mathrm{diag}]_{i,j,l,n} =
%>   \int_0^1 \hat{\varphi}_i \circ \hat{\mathbf{\gamma}}_n(s) 
%>   \hat{\varphi}_l \circ \hat{\mathbf{\gamma}}_n(s) 
%>   \hat{\varphi}_j\circ \hat{\mathbf{\gamma}}_n(s) \mathrm{d}s \,,
%> @f]
%> where the mapping @f$\hat{\mathbf{\gamma}}_n@f$ is defined in 
%> <code>gammaMap()</code>. The entries of matrix 
%> @f$\hat{\mathsf{{R}}}^\mathrm{offdiag} \in
%>    \mathbb{R}^{N\times N\times N\times 3\times 3}@f$
%> are given by
%> @f[
%> [\hat{\mathsf{{R}}}^\mathrm{offdiag}]_{i,j,l,n^-,n^+} =
%>   \int_0^1 \hat{\varphi}_i \circ \hat{\mathbf{\gamma}}_{n^-}(s) 
%>   \hat{\varphi}_l\circ \hat{\mathbf{\gamma}}_{n^-}(s)
%>   \hat{\varphi}_j\circ \hat{\mathbf{\vartheta}}_{n^-n^+} \circ
%>   \hat{\mathbf{\gamma}}_{n^-}(s) \mathrm{d}s \,,
%> @f]
%> with the mapping @f$\hat{\mathbf{\vartheta}}_{n^-n^+}@f$ as described in
%> <code>theta()</code>.
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0Tint <code>logical</code> arrays that mark each triangles
%>                    (interior) edges on which the matrix blocks should be
%>                    assembled @f$[K \times 3]@f$
%> @param refEdgePhiIntPhiIntPhiInt  Local matrix 
%>                    @f$\hat{\mathsf{{R}}}^\text{diag}@f$ as provided
%>                    by <code>integrateRefEdgePhiIntPhiIntPhiInt()</code>.
%>                    @f$[N \times N \times N_\mathrm{data} \times 3]@f$
%> @param refEdgePhiIntPhiExtPhiInt  Local matrix 
%>                    @f$\hat{\mathsf{{R}}}^\text{offdiag}@f$ as provided
%>                    by <code>integrateRefEdgePhiIntPhiIntPhiExt()</code>
%>                    using <code>permute()</code>.
%>                    @f$[N \times N \times N_\mathrm{data} \times 3 \times 3]@f$
%> @param dataDisc    A representation of the discrete function 
%>                    @f$d_h(\mathbf(x))@f$, e.g., as computed by 
%>                    <code>projectFuncCont2DataDisc()</code>
%>                    @f$[K \times N_\mathrm{data}]@f$
%> @param areaNuE0Tint (optional) argument to provide precomputed values
%>                    for the products of <code>markE0Tint</code>,
%>                    <code>g.areaE0T</code>, and <code>g.nuE0T</code>
%>                    @f$[3 \times 2 \text{ cell}]@f$
%> @retval ret        The assembled matrices @f$[2 \times 1 \text{ cell}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2015 Florian Frank, Balthasar Reuter, Vadym Aizinger
%> Modified by Hennes Hajduk, 2016-04-06
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
function ret = assembleMatEdgePhiPhiFuncDiscIntNu(g, markE0Tint, refEdgePhiIntPhiIntPhiInt, refEdgePhiIntPhiExtPhiInt, dataDisc, areaNuE0Tint)
[K, dataN] = size(dataDisc);
N = size(refEdgePhiIntPhiIntPhiInt, 1);

% Check function arguments that are directly used
validateattributes(dataDisc, {'numeric'}, {'size', [g.numT dataN]}, mfilename, 'dataDisc');
validateattributes(markE0Tint, {'logical'}, {'size', [K 3]}, mfilename, 'markE0Tint');
validateattributes(refEdgePhiIntPhiIntPhiInt, {'numeric'}, {'size', [N N dataN 3]}, mfilename, 'refEdgePhiIntPhiIntPhiInt');
validateattributes(refEdgePhiIntPhiExtPhiInt, {'numeric'}, {'size', [N N dataN 3 3]}, mfilename, 'refEdgePhiIntPhiExtPhiInt');

if nargin > 5
  if isfield(g, 'areaNuE0TE0T')
    ret = assembleMatEdgePhiPhiFuncDiscNu_withAreaNuE0Tint_wiAreaNuE0TE0T(g, refEdgePhiIntPhiIntPhiInt, refEdgePhiIntPhiExtPhiInt, dataDisc, areaNuE0Tint);
  else
    ret = assembleMatEdgePhiPhiFuncDiscNu_withAreaNuE0Tint_noAreaNuE0TE0T(g, refEdgePhiIntPhiIntPhiInt, refEdgePhiIntPhiExtPhiInt, dataDisc, areaNuE0Tint);
  end
else
  if isfield(g, 'areaNuE0TE0T')
    if isfield(g, 'areaNuE0T')
      ret = assembleMatEdgePhiPhiFuncDiscNu_withAreaNuE0TE0T_withAreaNuE0T(g, markE0Tint, refEdgePhiIntPhiIntPhiInt, refEdgePhiIntPhiExtPhiInt, dataDisc);
    else
      ret = assembleMatEdgePhiPhiFuncDiscNu_withAreaNuE0TE0T_noAreaNuE0T(g, markE0Tint, refEdgePhiIntPhiIntPhiInt, refEdgePhiIntPhiExtPhiInt, dataDisc);
    end
  else
    if isfield(g, 'areaNuE0T')
      ret = assembleMatEdgePhiPhiFuncDiscNu_noAreaNuE0TE0T_withAreaNuE0T(g, markE0Tint, refEdgePhiIntPhiIntPhiInt, refEdgePhiIntPhiExtPhiInt, dataDisc);
    else
      ret = assembleMatEdgePhiPhiFuncDiscNu_noAreaNuE0TE0T_noAreaNuE0T(g, markE0Tint, refEdgePhiIntPhiIntPhiInt, refEdgePhiIntPhiExtPhiInt, dataDisc);
    end    
  end
end % if
end % function
%
%===============================================================================
%> @brief Helper function for the case that assembleMatEdgePhiPhiFuncDiscNu()
%> was called with a precomputed field areaNuE0Tint and parameter g provides
%> a precomputed field areaNuE0TE0T.
%
function ret = assembleMatEdgePhiPhiFuncDiscNu_withAreaNuE0Tint_wiAreaNuE0TE0T(g, refEdgePhiIntPhiIntPhiInt, refEdgePhiIntPhiExtPhiInt, dataDisc, areaNuE0Tint)
[K, dataN] = size(dataDisc);
N = size(refEdgePhiIntPhiIntPhiInt, 1);

% Assemble matrices
ret = cell(2,1); 
ret{1} = sparse(K*N, K*N); 
ret{2} = sparse(K*N, K*N);
for nn = 1 : 3
  % Off-diagonal blocks
  for np = 1 : 3
    RtildeT = zeros(K*N, N);
    for l = 1 : dataN
      RtildeT = RtildeT + kron(dataDisc(:,l), refEdgePhiIntPhiExtPhiInt(:,:,l,nn,np));
    end % for
    ret{1} = ret{1} + 0.5 * kronVec(g.areaNuE0TE0T{nn,np,1}, RtildeT);
    ret{2} = ret{2} + 0.5 * kronVec(g.areaNuE0TE0T{nn,np,2}, RtildeT);
  end % for
  % Diagonal blocks
  for l = 1 : dataN
    ret{1} = ret{1} + kron(spdiags(0.5 * areaNuE0Tint{nn,1} .* dataDisc(:,l),0,K,K), refEdgePhiIntPhiIntPhiInt(:,:,l,nn));
    ret{2} = ret{2} + kron(spdiags(0.5 * areaNuE0Tint{nn,2} .* dataDisc(:,l),0,K,K), refEdgePhiIntPhiIntPhiInt(:,:,l,nn));
  end % for
end % for
end % function
%
%===============================================================================
%> @brief Helper function for the case that assembleMatEdgePhiPhiFuncDiscNu()
%> was called with a precomputed field areaNuE0Tint and parameter g provides
%> no precomputed field areaNuE0TE0T.
%
function ret = assembleMatEdgePhiPhiFuncDiscNu_withAreaNuE0Tint_noAreaNuE0TE0T(g, refEdgePhiIntPhiIntPhiInt, refEdgePhiIntPhiExtPhiInt, dataDisc, areaNuE0Tint)
[K, dataN] = size(dataDisc);
N = size(refEdgePhiIntPhiIntPhiInt, 1);

% Assemble matrices
ret = cell(2,1); 
ret{1} = sparse(K*N, K*N); 
ret{2} = sparse(K*N, K*N);
for nn = 1 : 3
  % Off-diagonal blocks
  Rkn = 0.5 * g.areaE0T(:,nn);
  for np = 1 : 3
    markE0TE0TtimesRkn1 = spdiags(Rkn .* g.nuE0T(:,nn,1), 0,K,K) * g.markE0TE0T{nn, np};
    markE0TE0TtimesRkn2 = spdiags(Rkn .* g.nuE0T(:,nn,2), 0,K,K) * g.markE0TE0T{nn, np};
    RtildeT = zeros(K*N, N);
    for l = 1 : dataN
      RtildeT = RtildeT + kron(dataDisc(:,l), refEdgePhiIntPhiExtPhiInt(:,:,l,nn,np));
    end % for
    ret{1} = ret{1} + kronVec(markE0TE0TtimesRkn1, RtildeT);
    ret{2} = ret{2} + kronVec(markE0TE0TtimesRkn2, RtildeT);
  end % for
  % Diagonal blocks
  for l = 1 : dataN
    ret{1} = ret{1} + kron(spdiags(0.5 * areaNuE0Tint{nn,1} .* dataDisc(:,l),0,K,K), refEdgePhiIntPhiIntPhiInt(:,:,l,nn));
    ret{2} = ret{2} + kron(spdiags(0.5 * areaNuE0Tint{nn,2} .* dataDisc(:,l),0,K,K), refEdgePhiIntPhiIntPhiInt(:,:,l,nn));
  end % for
end % for
end % function
%
%===============================================================================
%> @brief Helper function for the case that assembleMatEdgePhiPhiFuncDiscNu()
%> was called with no precomputed field areaNuE0Tint and parameter g provides
%> precomputed fields areaNuE0TE0T and areaNuE0T.
%
function ret = assembleMatEdgePhiPhiFuncDiscNu_withAreaNuE0TE0T_withAreaNuE0T(g, markE0Tint, refEdgePhiIntPhiIntPhiInt, refEdgePhiIntPhiExtPhiInt, dataDisc)
[K, dataN] = size(dataDisc);
N = size(refEdgePhiIntPhiIntPhiInt, 1);

% Assemble matrices
ret = cell(2,1); 
ret{1} = sparse(K*N, K*N); 
ret{2} = sparse(K*N, K*N);
for nn = 1 : 3
  % Off-diagonal blocks
  for np = 1 : 3
    RtildeT = zeros(K*N, N);
    for l = 1 : dataN
      RtildeT = RtildeT + kron(dataDisc(:,l), refEdgePhiIntPhiExtPhiInt(:,:,l,nn,np));
    end % for
    ret{1} = ret{1} + 0.5 * kronVec(g.areaNuE0TE0T{nn,np,1}, RtildeT);
    ret{2} = ret{2} + 0.5 * kronVec(g.areaNuE0TE0T{nn,np,2}, RtildeT);
  end % for
  % Diagonal blocks
  for l = 1 : dataN
    ret{1} = ret{1} + kron(spdiags(0.5 * g.areaNuE0T{nn,1} .* markE0Tint(:, nn) .* dataDisc(:,l), 0,K,K), refEdgePhiIntPhiIntPhiInt(:,:,l,nn));
    ret{2} = ret{2} + kron(spdiags(0.5 * g.areaNuE0T{nn,2} .* markE0Tint(:, nn) .* dataDisc(:,l), 0,K,K), refEdgePhiIntPhiIntPhiInt(:,:,l,nn));
  end % for
end % for
end % function
%
%===============================================================================
%> @brief Helper function for the case that assembleMatEdgePhiPhiFuncDiscNu()
%> was called with no precomputed field areaNuE0Tint and parameter g provides
%> a precomputed field areaNuE0TE0T but no field areaNuE0T.
%
function ret = assembleMatEdgePhiPhiFuncDiscNu_withAreaNuE0TE0T_noAreaNuE0T(g, markE0Tint, refEdgePhiIntPhiIntPhiInt, refEdgePhiIntPhiExtPhiInt, dataDisc)
[K, dataN] = size(dataDisc);
N = size(refEdgePhiIntPhiIntPhiInt, 1);

% Assemble matrices
ret = cell(2,1); 
ret{1} = sparse(K*N, K*N); 
ret{2} = sparse(K*N, K*N);
for nn = 1 : 3
  % Off-diagonal blocks
  for np = 1 : 3
    RtildeT = zeros(K*N, N);
    for l = 1 : dataN
      RtildeT = RtildeT + kron(dataDisc(:,l), refEdgePhiIntPhiExtPhiInt(:,:,l,nn,np));
    end % for
    ret{1} = ret{1} + 0.5 * kronVec(g.areaNuE0TE0T{nn,np,1}, RtildeT);
    ret{2} = ret{2} + 0.5 * kronVec(g.areaNuE0TE0T{nn,np,2}, RtildeT);
  end % for
  % Diagonal blocks
  Rkn = 0.5 * g.areaE0T(:,nn) .* markE0Tint(:, nn);
  for l = 1 : dataN
    ret{1} = ret{1} + kron(spdiags(Rkn .* g.nuE0T(:,nn,1) .* dataDisc(:,l), 0,K,K), refEdgePhiIntPhiIntPhiInt(:,:,l,nn));
    ret{2} = ret{2} + kron(spdiags(Rkn .* g.nuE0T(:,nn,2) .* dataDisc(:,l), 0,K,K), refEdgePhiIntPhiIntPhiInt(:,:,l,nn));
  end % for
end % for
end % function
%
%===============================================================================
%> @brief Helper function for the case that assembleMatEdgePhiPhiFuncDiscNu()
%> was called with no precomputed field areaNuE0Tint and parameter g provides
%> a precomputed field areaNuE0T but no field areaNuE0TE0T.
%
function ret = assembleMatEdgePhiPhiFuncDiscNu_noAreaNuE0TE0T_withAreaNuE0T(g, markE0Tint, refEdgePhiIntPhiIntPhiInt, refEdgePhiIntPhiExtPhiInt, dataDisc)
[K, dataN] = size(dataDisc);
N = size(refEdgePhiIntPhiIntPhiInt, 1);

% Assemble matrices
ret = cell(2,1); 
ret{1} = sparse(K*N, K*N); 
ret{2} = sparse(K*N, K*N);
for nn = 1 : 3
  % Off-diagonal blocks
  Rkn = 0.5 * g.areaE0T(:,nn);
  for np = 1 : 3
    markE0TE0TtimesRkn1 = spdiags(Rkn .* g.nuE0T(:,nn,1), 0,K,K) * g.markE0TE0T{nn, np};
    markE0TE0TtimesRkn2 = spdiags(Rkn .* g.nuE0T(:,nn,2), 0,K,K) * g.markE0TE0T{nn, np};
    RtildeT = zeros(K*N, N);
    for l = 1 : dataN
      RtildeT = RtildeT + kron(dataDisc(:,l), refEdgePhiIntPhiExtPhiInt(:,:,l,nn,np));
    end % for
    ret{1} = ret{1} + kronVec(markE0TE0TtimesRkn1, RtildeT);
    ret{2} = ret{2} + kronVec(markE0TE0TtimesRkn2, RtildeT);
  end % for
  % Diagonal blocks
  for l = 1 : dataN
    ret{1} = ret{1} + kron(spdiags(0.5 * g.areaNuE0T{nn,1} .* markE0Tint(:, nn) .* dataDisc(:,l), 0,K,K), refEdgePhiIntPhiIntPhiInt(:,:,l,nn));
    ret{2} = ret{2} + kron(spdiags(0.5 * g.areaNuE0T{nn,2} .* markE0Tint(:, nn) .* dataDisc(:,l), 0,K,K), refEdgePhiIntPhiIntPhiInt(:,:,l,nn));
  end % for
end % for
end % function
%
%===============================================================================
%> @brief Helper function for the case that assembleMatEdgePhiPhiFuncDiscNu()
%> was called with no precomputed field areaNuE0Tint and parameter g provides
%> no precomputed fields areaNuE0TE0T or areaNuE0T.
%
function ret = assembleMatEdgePhiPhiFuncDiscNu_noAreaNuE0TE0T_noAreaNuE0T(g, markE0Tint, refEdgePhiIntPhiIntPhiInt, refEdgePhiIntPhiExtPhiInt, dataDisc)
[K, dataN] = size(dataDisc);
N = size(refEdgePhiIntPhiIntPhiInt, 1);

% Assemble matrices
ret = cell(2, 1); ret{1} = sparse(K*N, K*N); ret{2} = sparse(K*N, K*N);
for nn = 1 : 3
  % Off-diagonal blocks
  Rkn = 0.5 * g.areaE0T(:,nn);
  for np = 1 : 3
    markE0TE0TtimesRkn1 = spdiags(Rkn .* g.nuE0T(:,nn,1), 0,K,K) * g.markE0TE0T{nn, np};
    markE0TE0TtimesRkn2 = spdiags(Rkn .* g.nuE0T(:,nn,2), 0,K,K) * g.markE0TE0T{nn, np};
    RtildeT = zeros(K*N, N);
    for l = 1 : dataN
      RtildeT = RtildeT + kron(dataDisc(:,l), refEdgePhiIntPhiExtPhiInt(:,:,l,nn,np));
    end % for
    ret{1} = ret{1} + kronVec(markE0TE0TtimesRkn1, RtildeT);
    ret{2} = ret{2} + kronVec(markE0TE0TtimesRkn2, RtildeT);
  end % for
  % Diagonal blocks
  Rkn = Rkn .* markE0Tint(:,nn);
  for l = 1 : dataN
    ret{1} = ret{1} + kron(spdiags(Rkn .* g.nuE0T(:,nn,1) .* dataDisc(:,l), 0,K,K), refEdgePhiIntPhiIntPhiInt(:,:,l,nn));
    ret{2} = ret{2} + kron(spdiags(Rkn .* g.nuE0T(:,nn,2) .* dataDisc(:,l), 0,K,K), refEdgePhiIntPhiIntPhiInt(:,:,l,nn));
  end % for
end % for
end % function