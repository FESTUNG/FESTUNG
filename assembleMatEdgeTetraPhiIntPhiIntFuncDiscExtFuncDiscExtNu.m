% Assembles two matrices containing integrals over edges of products of
% two basis function with two discontinuous coefficient functions
% from the neighboring element and with a component of the edge normal.

%===============================================================================
%> @file assembleMatEdgeTetraPhiIntPhiIntFuncDiscExtFuncDiscExtNu.m
%>
%> @brief Assembles two matrices containing integrals over edges of products of
%>        two basis function with two discontinuous coefficient functions
%>        from the neighboring element and with a component of the edge normal.
%===============================================================================
%>
%> @brief Assembles two matrices containing integrals over edges of products of
%>        two basis function with two discontinuous coefficient functions
%>        from the neighboring element and with a component of the edge normal.
%>
%> The matrices @f$\mathsf{{R}}^m\in\mathbb{R}^{KN\times KN}@f$ are defined as
%> @f[
%> [\mathsf{{R}}^m]_{(k^--1)N+i,(k^+-1)N+j} =
%>   \nu^m_{k^-n^-} \sum_{l=1}^N D_{k^+l}(t) \sum_{ll=1}^N F_{k^+ll}(t) 
%>     \int_{E_{k^-n^-}} \varphi_{k^-i} \varphi_{k^-j} \varphi_{k^+l} 
%>     \varphi_{k^+ll} \mathrm{d}s \,.
%> @f]
%> with @f$\nu_{kn}@f$ the @f$m@f$-th component (@f$m\in\{1,2\}@f$) of the unit
%> normal and @f$D_{kl}(t)@f$ and @f$F_{kll}(t)@f$ the coefficients of the 
%> discrete representation of functions
%> @f$ d_h(t, \mathbf{x}) = \sum_{l=1}^N D_{kl}(t) \varphi_{kl}(\mathbf{x}) @f$
%> and
%> @f$ f_h(t, \mathbf{x}) = \sum_{l=1}^N F_{kl}(t) \varphi_{kl}(\mathbf{x}) @f$
%> on an element @f$T_k@f$.
%> All other entries are zero.
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>domainRectTrap()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0T    <code>logical</code> arrays that mark each elements
%>                    edges on which the matrix blocks should be
%>                    assembled @f$[K \times 4]@f$
%> @param refEdgePhiIntPhiIntPhiExtPhiExt  Local matrix 
%>                    @f$\hat{\mathsf{{R}}}^\text{offdiag}@f$ as provided
%>                    by <code>integrateRefEdgeTetraPhiIntPhiIntPhiExtPhiExt()</code>.
%>                    @f$[N \times N \times N \times N \times 4]@f$
%> @param dataDisc1    A representation of the discrete function 
%>                    @f$d_h(\mathbf(x))@f$, e.g., as computed by 
%>                    <code>projectFuncCont2DataDiscTetra()</code>
%>                    @f$[K \times N]@f$
%> @param dataDisc2    A representation of the discrete function 
%>                    @f$f_h(\mathbf(x))@f$, e.g., as computed by 
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
function ret = assembleMatEdgeTetraPhiIntPhiIntFuncDiscExtFuncDiscExtNu(g, markE0T, refEdgePhiIntPhiIntPhiExtPhiExt, dataDisc1, dataDisc2)
[K, N] = size(dataDisc1);
ret = { sparse(K*N, K*N), sparse(K*N, K*N) };
for n = 1 : 4
  RtildeT1 = zeros(K*N, N);
  for l = 1 : N
    RtildeT2 = zeros(K*N, N);
    for ll = 1 : N
      RtildeT2 = RtildeT2 + kron(dataDisc2(:,ll), refEdgePhiIntPhiIntPhiExtPhiExt(:,:,l,ll,n));
    end % for ll
    RtildeT1 = RtildeT1 + kron(dataDisc1(:,l), ones(N, N)) .* RtildeT2;
  end % for l
  for m = 1 : 2
    markAreaNuE0TE0T = spdiags(g.areaE0T(:,n) .* g.nuE0T(:,n,m) .* markE0T(:,n), 0, K, K) * g.markE0TE0T{n};
    ret{m} = ret{m} + kronVec(markAreaNuE0TE0T.', RtildeT1).';
  end % for m
end  % for n
end % function
