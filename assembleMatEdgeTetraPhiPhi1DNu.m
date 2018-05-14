% Assembles two matrices containing integrals over interior edges of products of
% a two-dimensional and a one-dimensional basis function with a component of the
% edge normal.

%===============================================================================
%> @file
%>
%> @brief Assembles two matrices containing integrals over interior edges of 
%>        products of a two-dimensional and a one-dimensional basis function 
%>        with a component of the edge normal.
%===============================================================================
%>
%> @brief Assembles two matrices containing integrals over interior edges of 
%>        products of a two-dimensional and a one-dimensional basis function 
%>        with a component of the edge normal.
%>
%> The matrices @f$\mathsf{{Q}}^m \in 
%>      \mathbb{R}^{KN\times \overline{K}\overline{N}}@f$
%> consist of two kinds of contributions: diagonal blocks and off-diagonal 
%> blocks. Diagonal blocks are defined as 
%> @f[
%> [\mathsf{{Q}}^m]_{(k-1)N+i,(\overline{k}-1)\overline{N}+j} = 
%>  \frac{1}{2} \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_\Omega}
%>  \nu_{kn}^m  \int_{E_{kn}} \varphi_{ki} \phi_{\overline{k}j} \mathrm{d}s \,.
%> @f]
%> and off-diagonal blocks are defined as
%> @f[
%> [\mathsf{{Q}}^m]_{(k^--1)N+i,(\overline{k}^+-1)\overline{N}+j} = 
%>  \frac{1}{2} \nu_{k^-n^-}^m \int_{E_{k^-n^-}} 
%>   \varphi_{k^-i} \phi_{\overline{k}^+j} \mathrm{d}s \,.
%> @f]
%> with @f$\nu_{kn}@f$ the @f$m@f$-th component (@f$m\in\{1,2\}@f$) of the edge
%> normal.
%> Entries in off-diagonal blocks are only non-zero for pairs of elements
%> @f$T_{k^-}, T_{k^+}@f$ with @f$\partial T_{k^-} \cap T_{k^+} \ne\emptyset@f$.
%> Note that the local edge index @f$n^-@f$ is given implicitly, since 
%> @f$\partial T_{k^-} \cap T_{k^+} \ne\emptyset@f$ consist of exactly one
%> edge @f$E_{k^-n^-} = E_{k^+n^+}@f$.
%> To allow for vectorization, the assembly is reformulated as
%> @f$\mathsf{Q}^m = \mathsf{Q}^{m,\mathrm{diag}} + 
%>    \mathsf{Q}^{m,\mathrm{offdiag}}@f$ with the blocks defined as
%> @f[
%> \mathsf{{Q}}^{m,\mathrm{diag}} = \frac{1}{2} \sum_{n=1}^4
%>   \begin{bmatrix}
%>     \delta_{E_{1n}\in\mathcal{E}_\Omega} &   & \\
%>     & ~\ddots~ & \\
%>     &          & \delta_{E_{Kn}\in\mathcal{E}_\Omega}
%>   \end{bmatrix} \circ \begin{bmatrix}
%>     \nu^m_{1n} | E_{1n} | &   & \\
%>     & ~\ddots~ & \\
%>     &          & \nu^m_{Kn} | E_{Kn} |
%>   \end{bmatrix} \otimes [\hat{\mathsf{Q}}^\mathrm{diag}]_{:,:,n}\;,
%> @f]
%> and
%> @f[
%> \mathsf{Q}^{m,\mathrm{offdiag}} = \frac{1}{2}\sum_{n^-=1}^4
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
%> [\hat{\mathsf{Q}}^\mathrm{offdiag}]_{:,:,n^-}\;,
%> @f]
%> where @f$\delta_{E_{kn}\in\mathcal{E}_\Omega}, \delta_{E_{in^-}=E_{jn^+}}@f$ 
%> denote the Kronecker delta, @f$\circ@f$ denotes the Hadamard product, and 
%> @f$\otimes@f$ denotes the Kronecker product.
%> The index @f$n^+@f$ is given implicitly by @f$n^-@f$ as described in
%> <code>mapLocalEdgeTetra()</code>.
%>
%> The entries of matrix 
%> @f$\hat{\mathsf{Q}}^\mathrm{diag}\in\mathbb{R}^{N\times \overline{N}\times4}@f$
%> are given by
%> @f[
%> [\hat{\mathsf{Q}}^\mathrm{diag}]_{i,j,n} =
%>   \int_0^1 \hat{\varphi}_i \circ \hat{\mathbf{\gamma}}_n(s) 
%>   \hat{\phi}_j\circ [\hat{\mathbf{\gamma}}_n(s)]_1 \mathrm{d}s \,,
%> @f]
%> where the mapping @f$\hat{\mathbf{\gamma}}_n@f$ is defined in 
%> <code>gammaMapTetra()</code>. The entries of matrix
%> @f$\hat{\mathsf{Q}}^\mathrm{offdiag} \in 
%>    \mathbb{R}^{N\times \overline{N}\times 4}@f$ are defined as
%> @f[
%> [\hat{\mathsf{Q}}^\mathrm{offdiag}]_{i,j,n^-,n^+} =
%>   \int_0^1 \hat{\varphi}_i \circ \hat{\mathbf{\gamma}}_{n^-}(s) 
%>   \hat{\phi}_j\circ [\hat{\mathbf{\gamma}}_{n^+}(s)]_1 \mathrm{d}s \,.
%> @f]
%>
%> @param  g2D        The lists describing the geometric and topological 
%>                    properties of a quadrilateral triangulation (see 
%>                    <code>domainRectTrap()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  g1D        The lists describing the geometric and topological 
%>                    properties of a matching 1D triangulation (see 
%>                    <code>generateGridData1D()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0T    <code>logical</code> arrays that mark each elements
%>                    edges on which the matrix blocks should be
%>                    assembled @f$[K \times 4]@f$
%> @param refEdgePhiIntPhi1DInt  Local matrix 
%>                    @f$\hat{\mathsf{Q}}^\text{diag}@f$ as provided
%>                    by <code>integrateRefEdgeTetraPhiIntPhi1DInt()</code>.
%>                    @f$[N \times N \times 4]@f$
%> @param refEdgePhiIntPhi1DExt Local matrix 
%>                    @f$\hat{\mathsf{Q}}^\text{offdiag}@f$ as provided
%>                    by <code>integrateRefEdgeTetraPhiIntPhi1DExt()</code>.
%>                    @f$[N \times N \times 4]@f$
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
function ret = assembleMatEdgeTetraPhiPhi1DNu(g2D, g1D, markE0T, refEdgePhiIntPhi1DInt, refEdgePhiIntPhi1DExt)
K = g2D.numT; barK = g1D.numT;
[N, barN, ~] = size(refEdgePhiIntPhi1DInt);
ret = { sparse(K*N, barK*barN), sparse(K*N, barK*barN) };
for n = 1 : 4
  areaE0T = 0.5 * markE0T(:,n) .* g2D.areaE0T(:,n);
  for m = 1 : 2
    areaNuE0T = areaE0T .* g2D.nuE0T(:,n,m);
    ret{m} = ret{m} + ...
             kron(bsxfun(@times, g1D.markT2DT, areaNuE0T), refEdgePhiIntPhi1DInt(:,:,n)) + ...
             kron(bsxfun(@times, g2D.markE0TE0T{n} * double(g1D.markT2DT), areaNuE0T), refEdgePhiIntPhi1DExt(:,:,n));
  end % for m
end  % for n
end  % function