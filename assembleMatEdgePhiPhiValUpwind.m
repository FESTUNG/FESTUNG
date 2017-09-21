% Assembles a matrix containing integrals over interior edges of products of
% two basis functions with the upwind value of a, for each quadrature point 
% specified, function.

%===============================================================================
%> @file assembleMatEdgePhiPhiValUpwind.m
%>
%> @brief Assembles a matrix containing integrals over interior edges of products of
%>        two basis functions with the upwind value of a, for each quadrature 
%>        point specified, function.
%===============================================================================
%>
%> @brief Assembles matrix @f$\mathsf{R}@f$ containing integrals over
%>        edges of products of two basis functions with the upwind value of 
%>        a, for each quadrature point specified, function.
%>
%> The matrix @f$\mathsf{R} \in \mathbb{R}^{KN\times KN}@f$
%> consist of two kinds of contributions: diagonal blocks and off-diagonal 
%> blocks. Diagonal blocks are defined as 
%> @f[
%> [\mathsf{R}]_{(k-1)N+i,(k-1)N+j} = 
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_\Omega}
%>  \int_{E_{kn}} \varphi_{ki} \varphi_{kj} 
%>  d_h(\mathbf{x}) \delta_{d_h(\mathbf{x}) \ge 0} \mathrm{d}s \,.
%> @f]
%> and off-diagonal blocks are defined as
%> @f[
%> [\mathsf{R}]_{(k^--1)N+i,(k^+-1)N+j} = 
%>  \int_{E_{k^-n^-}} \varphi_{k^-i} \varphi_{k^+j}
%>  d_h(\mathbf{x}) \delta_{d_h(\mathbf{x}) \ge 0} \mathrm{d}s \,.
%> @f]
%> with @f$d_h(\mathbf{x})@f$ a coefficient function and
%> @f$\delta_{d_h(\mathbf{x}) \ge 0}@f$ the Kronecker delta.
%> Entries in off-diagonal blocks are only non-zero for pairs of triangles
%> @f$T_{k^-}, T_{k^+}@f$ with @f$\partial T_{k^-} \cap T_{k^+} \ne\emptyset@f$.
%> Note that the local edge index @f$n^-@f$ is given implicitly, since 
%> @f$\partial T_{k^-} \cap T_{k^+} \ne\emptyset@f$ consist of exactly one
%> edge @f$E_{k^-n^-} = E_{k^+n^+}@f$.
%>
%> The coefficient function @f$d_h(\mathbf{x})@f$ can, e.g., be the
%> normal velocity @f$\mathbf{u}(\mathbf{x})\cdot\mathbf{\nu}_{k^-n^-}@f$, 
%> thus allowing a vectorized evaluation of an upwind flux formulation.
%> 
%> To allow for vectorization, the assembly is reformulated as
%> @f$\mathsf{R} = \mathsf{R}^\mathrm{diag} + 
%>    \mathsf{R}^\mathrm{offdiag}@f$ with the blocks defined as
%> @f[
%> \mathsf{R}^\mathrm{diag} = \sum_{n=1}^3 \sum_{r=1}^R \omega_r
%>   \begin{bmatrix}
%>     \delta_{E_{1n}\in\mathcal{E}_\Omega} &   & \\
%>     & ~\ddots~ & \\
%>     &          & \delta_{E_{Kn}\in\mathcal{E}_\Omega}
%>   \end{bmatrix} \circ \begin{bmatrix}
%>     | E_{1n} | &   & \\
%>     & ~\ddots~ & \\
%>     &          & | E_{Kn} |
%>   \end{bmatrix} \circ \begin{bmatrix}
%>     d_h \circ \mathbf{F}_1 \circ \hat{\gamma}_n(\hat{q}_r) \delta_{d_h \ge 0} &   & \\
%>     & ~\ddots~ & \\
%>     &          & d_h \circ \mathbf{F}_K \circ \hat{\gamma}_n(\hat{q}_r) \delta_{d_h \ge 0}
%>   \end{bmatrix} \otimes [\hat{\mathsf{R}}^\mathrm{diag}]_{:,:,n,r}\;,
%> @f]
%> and
%> @f[
%> \mathsf{R}^\mathrm{offdiag} = \sum_{n^-=1}^3\sum_{n^+=1}^3 \sum_{r=1}^R \omega_r
%>   \begin{bmatrix}
%>     0&\delta_{E_{1n^-} = E_{2n^+}}&\dots&\dots&\delta_{E_{1n^-}=E_{Kn^+}} \\
%>     \delta_{E_{2n^-} = E_{1n^+}}&0&\ddots& &\vdots \\
%>     \vdots & \ddots & \ddots & \ddots & \vdots \\
%>     \vdots & & \ddots & 0 & \delta_{E_{(K-1)n^-}=E_{Kn^+}} \\
%>     \delta_{E_{Kn^-} = E_{1n^+}}&\dots&\dots&\delta_{E_{Kn^-} = E_{(K-1)n^+}} &0
%>   \end{bmatrix} \circ \begin{bmatrix}
%>     d_h \circ \mathbf{F}_1 \circ \hat{\gamma}_{n^-}(\hat{q}_r)\delta_{d_h\ge0} |E_{1n^-}| & \dots & d_h \circ \mathbf{F}_1 \circ \hat{\gamma}_{n^-}(\hat{q}_r)\delta_{d_h\ge0} |E_{1n^-}| \\
%>     \vdots & & \vdots \\
%>     d_h \circ \mathbf{F}_K \circ \hat{\gamma}_{n^-}(\hat{q}_r)\delta_{d_h\ge0} |E_{Kn^-}|  & \dots & d_h \circ \mathbf{F}_K \circ \hat{\gamma}_{n^-}(\hat{q}_r)\delta_{d_h\ge0} |E_{Kn^-}|  
%>   \end{bmatrix} \otimes 
%> [\hat{\mathsf{R}}^\mathrm{offdiag}]_{:,:,n^-,n^+,r}\;,
%> @f]
%> where @f$\delta_{E_{kn}\in\mathcal{E}_\Omega}, \delta_{E_{in^-}=E_{jn^+}}@f$ 
%> denote the Kronecker delta, @f$\circ@f$ denotes the Hadamard product, 
%> @f$\otimes@f$ denotes the Kronecker product, 
%> @f$\mathbf{F}_k: \hat{T} \ni \hat{\mathbf{x}} \rightarrow \mathbf{x} \in T_k@f$.
%> is the affine mapping from reference triangle @f$\hat{T}@f$ to the
%> physical triangle @f$T_k@f$, and
%> @f$\hat{\mathbf{\gamma}}_n: [0,1] \ni s \rightarrow \hat{\mathbf{x}} \in \hat{T}@f$
%> maps from the refence interval to the nth edge in the reference triangle and is
%> defined in <code>gammaMap()</code>.
%> Quadrature weights @f$\omega_r@f$ and points @f$\hat{q}_r@f$ are obtained
%> from <code>quadRule1D()</code>.
%>
%> The entries of matrix 
%> @f$\hat{\mathsf{R}}^\mathrm{diag}\in\mathbb{R}^{N\times N\times 3\times R}@f$
%> are given by
%> @f[
%> [\hat{\mathsf{R}}^\mathrm{diag}]_{i,j,n,r} =
%>   \hat{\varphi}_i \circ \hat{\mathbf{\gamma}}_n(\hat{q}_r) 
%>   \hat{\varphi}_j\circ \hat{\mathbf{\gamma}}_n(\hat{q}_r) \,.
%> @f]
%> The entries of matrix
%> @f$\hat{\mathsf{R}}^\mathrm{offdiag} \in 
%>    \mathbb{R}^{N\times N\times 3\times 3\times R}@f$ are defined as
%> @f[
%> [\hat{\mathsf{R}}^\mathrm{offdiag}]_{i,j,n^-,n^+,r} =
%>   \hat{\varphi}_i \circ \hat{\mathbf{\gamma}}_{n^-}(\hat{q}_r) 
%>   \hat{\varphi}_j\circ \hat{\mathbf{\vartheta}}_{n^-n^+} \circ
%>   \hat{\mathbf{\gamma}}_{n^-}(\hat{q}_r)\,,
%> @f]
%> with the mapping @f$\hat{\mathbf{\vartheta}}_{n^-n^+}@f$ as described in
%> <code>theta()</code>.
%>
%> To speed up computation, the off-diagonal blocks are reformulated using
%> the operation @f$\otimes_\mathrm{V}@f$ defined by <code>kronVec()</code>, and read now as
%> @f[
%> \mathsf{R}^\mathrm{offdiag} = \sum_{n^-=1}^3\sum_{n^+=1}^3
%> \begin{bmatrix}
%>  0                          & \delta_{E_{1n^-}=E_{2n^+}} & \hdots                & \hdots                & \delta_{E_{1n^-}=E_{Kn^+}} \\
%>  \delta_{E_{2n^-}=E_{1n^+}} & 0                          &   \ddots              &                       & \textstyle\vdots  \\ 
%>  \vdots                     &         \ddots             & \ddots                &          \ddots       & \vdots   \\
%>  \vdots                     & {}                         &    \ddots             & 0                     & \delta_{E_{(K-1)n^-}=E_{Kn^+}} \\
%>  \delta_{E_{Kn^-}=E_{1n^+}} &  \hdots                    & \hdots                & \delta_{E_{Kn^-}=E_{(K-1)n^+}}  & 0
%>  \end{bmatrix}
%>  \otimes_\mathrm{V} \left[ \tilde{\mathsf{R}}^\mathrm{offdiag} \right]_{:,:,n-,n+}
%> @f]
%> with
%> @f[
%> \left[ \tilde{\mathsf{R}}^\mathrm{offdiag} \right]_{(k-1)N+1:kN,:,n-,n+} =\;
%>  \sum_{r=1}^R \omega_r 
%>  \begin{bmatrix}
%>    |E_{1n^-}|d_h \circ \mathbf{F}_1 \circ \hat{\gamma}_{n^-}(\hat{q}_r)\delta_{d_h\ge0} \\
%>    \vdots \\
%>    |E_{Kn^-}|d_h \circ \mathbf{F}_K \circ \hat{\gamma}_{n^-}(\hat{q}_r)\delta_{d_h\ge0}
%>  \end{bmatrix}
%>  \otimes [\hat{\mathsf{R}}^\mathrm{offdiag}]_{:,:,n^-,n^+,r}\;.
%> @f]
%>
%> Note that this assembly routine can also includes terms from edges that are
%> part of an outflow boundary, i.e., 
%> @f$d_h(\mathbf{x} \ge 0, \mathbf{x}\in\partial\Omega@f$. 
%> Contributions from an inflow boundary must be treated separately, e.g.,
%> by <code>assembleVecEdgePhiIntFuncContVal()</code>.
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0Tbdr <code>logical</code> arrays that mark each triangles
%>                    (boundary) edges on which the matrix blocks should be
%>                    assembled @f$[K \times 3]@f$
%> @param refEdgePhiIntPhiIntOnQuad  Local matrix 
%>                    @f$\hat{\mathsf{R}}^\text{diag}@f$ as provided
%>                    by <code>integrateRefEdgePhiIntPhiIntPerQuad()</code>.
%>                    @f$[N \times N \times 3 \times R]@f$
%> @param refEdgePhiIntPhiExtOnQuad Local matrix 
%>                    @f$\hat{\mathsf{R}}^\text{offdiag}@f$ as provided
%>                    by <code>integrateRefEdgePhiIntPhiExtPerQuad()</code>.
%>                    @f$[N \times N \times 3 \times 3 \times R]@f$
%> @param valOnQuad   The function values 
%>                    @f$d_h \circ \mathbf{F}_k \circ \hat{\gamma}_n(\hat{q}_r)@f$.
%>                    For a normal velocity these can be computed by 
%>                    <code>computeFuncContNuOnQuadEdge()</code>.
%>                    @f$[K \times 3 \times R]@f$
%> @param areaE0Tbdr (optional) argument to provide precomputed values
%>                    for the products of <code>markE0Tbdr</code>, and
%>                    <code>g.areaE0T</code>
%>                    @f$[3 \times 1 \text{ cell}]@f$
%> @param  elem       (optional) <code>logical</code> arrays to provide the 
%>                    elements of the grid for which the computation is done.
%>                    @f$[K \times 3]@f$
%> @retval ret        The assembled matrix @f$[KN \times KN]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Florian Frank, Balthasar Reuter, Vadym Aizinger
%>                      Modified 09/02/16 by Hennes Hajduk
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
function ret = assembleMatEdgePhiPhiValUpwind(g, markE0Tbdr, refEdgePhiIntPhiIntOnQuad, refEdgePhiIntPhiExtOnQuad, valOnQuad, areaE0Tbdr, elem)
% Extract dimensions
[K, ~, R] = size(valOnQuad);
N = size(refEdgePhiIntPhiIntOnQuad, 1);

if nargin < 7
  elem = true(K,1);
end % if

% Check function arguments that are directly used
validateattributes(markE0Tbdr, {'logical'}, {'size', [K 3]}, mfilename, 'markE0Tbdr');
validateattributes(valOnQuad, {'numeric'}, {'size', [K 3 NaN]});
validateattributes(refEdgePhiIntPhiIntOnQuad, {'numeric'}, {'size', [N N 3 R]});
validateattributes(refEdgePhiIntPhiExtOnQuad, {'numeric'}, {'size', [N N 3 3 R]});
validateattributes(elem, {'logical'}, {'size', [K 1]}, mfilename, 'elem');

if nargin > 5 && ~isequal(areaE0Tbdr, [])
  ret = assembleMatEdgePhiPhiValUpwind_withAreaE0Tbdr(g, refEdgePhiIntPhiIntOnQuad, refEdgePhiIntPhiExtOnQuad, valOnQuad(elem,:,:), areaE0Tbdr, elem);
else
  ret = assembleMatEdgePhiPhiValUpwind_noAreaE0Tbdr(g, markE0Tbdr, refEdgePhiIntPhiIntOnQuad, refEdgePhiIntPhiExtOnQuad, valOnQuad(elem,:,:), elem);
end % for
end % function

function ret = assembleMatEdgePhiPhiValUpwind_withAreaE0Tbdr(g, refEdgePhiIntPhiIntOnQuad, refEdgePhiIntPhiExtOnQuad, valOnQuad, areaE0Tbdr, elem)
[K, ~, R] = size(valOnQuad);
N = size(refEdgePhiIntPhiIntOnQuad, 1);

% Assemble matrices
ret = sparse(K*N, g.numT*N);
retDiag = sparse(K*N, K*N);
for nn = 1 : 3
  % Diagonal blocks
  for r = 1 : R
    retDiag = retDiag + kron(spdiags(areaE0Tbdr{nn}(elem) .* valOnQuad(:, nn, r) .* (valOnQuad(:, nn, r) > 0), 0, K, K), refEdgePhiIntPhiIntOnQuad(:, :, nn, r));
  end % for
  % Off-diagonal blocks
  for np = 1 : 3
    RknTimesVal = sparse(K*N, N);
    for r = 1 : R
      RknTimesVal = RknTimesVal + kron(areaE0Tbdr{nn}(elem) .* valOnQuad(:, nn, r) .* sparse(valOnQuad(:, nn, r) < 0), refEdgePhiIntPhiExtOnQuad(:, :, nn, np, r));
    end % for
    ret = ret + kronVec(g.markE0TE0T{nn, np}(elem,:), RknTimesVal);
  end % for
end % for
ret(:,logical(kron(elem,true(N,1)))) = ret(:,logical(kron(elem,true(N,1)))) + retDiag;
end % function

function ret = assembleMatEdgePhiPhiValUpwind_noAreaE0Tbdr(g, markE0Tbdr, refEdgePhiIntPhiIntOnQuad, refEdgePhiIntPhiExtOnQuad, valOnQuad, elem)
[K, ~, R] = size(valOnQuad);
N = size(refEdgePhiIntPhiIntOnQuad, 1);

% Assemble matrices
ret = sparse(K*N, g.numT*N);
retDiag = sparse(K*N, K*N);
for nn = 1 : 3
  Rkn = markE0Tbdr(elem,nn) .* g.areaE0T(elem, nn);
  % Diagonal blocks
  for r = 1 : R
    retDiag = retDiag + kron(spdiags(Rkn .* valOnQuad(:, nn, r) .* (valOnQuad(:, nn, r) > 0), 0, K, K), refEdgePhiIntPhiIntOnQuad(:, :, nn, r));
  end % for
  % Off-diagonal blocks
  for np = 1 : 3
    RknTimesVal = sparse(K*N, N);
    for r = 1 : R
      RknTimesVal = RknTimesVal + kron(Rkn .* valOnQuad(:, nn, r) .* sparse(valOnQuad(:, nn, r) < 0), refEdgePhiIntPhiExtOnQuad(:, :, nn, np, r));
    end % for
    ret = ret + kronVec(g.markE0TE0T{nn, np}(elem,:), RknTimesVal);
  end % for
end % for
ret(:,logical(kron(elem,true(N,1)))) = ret(:,logical(kron(elem,true(N,1)))) + retDiag;
end % function
