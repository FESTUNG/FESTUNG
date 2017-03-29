% Assembles two matrices containing integrals over edges of products of
% two basis functions with a continuous coefficient function and with a
% component of the edge normal.

%===============================================================================
%> @file assembleMatEdgeTetraPhiIntPhiIntFuncContNu.m
%>
%> @brief Assembles two matrices containing integrals over edges of 
%>        products of two basis functions with a continuous coefficient 
%>        function and with a component of the edge normal.
%===============================================================================
%>
%> @brief Assembles the matrices @f$\mathsf{{R}}^m, m\in\{1,2\}@f$ containing 
%>        integrals over edges of products of two basis functions with 
%>        a continuous coefficient function and a component of the edge normal.
%>
%> The matrices @f$\mathsf{{R}}^{m,\mathrm{diag}} 
%>                  \in\mathbb{R}^{KN\times KN}@f$ are defined as
%> @f[
%> [\mathsf{{R}}^{m,\mathrm{diag}}]_{(k-1)N+i,(k-1)N+j} =
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_\Omega}
%>  \nu_{kn}^m \int_{E_kn} \varphi_{ki}\varphi_{kj}d(t) \,,
%> @f]
%> with @f$\nu_{kn}@f$ the @f$m@f$-th component (@f$m\in\{1,2\}@f$) of the unit
%> normal and @f$d(t)@f$ the continuous function.
%> All other entries are zero.
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>domainRectTrap()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0T    <code>logical</code> arrays that mark each elements
%>                    edges on which the matrix blocks should be
%>                    assembled @f$[K \times 4]@f$
%> @param funcCont    A function_handle for the continuous function.
%> @param qOrd        The order of the quadrature rule.
%> @param refEdgePhiIntPhiIntPerQuad  Local matrix 
%>                    @f$\hat{\mathsf{{R}}}^\text{diag}@f$ as provided
%>                    by <code>integrateRefEdgeTetraPhiIntPhiIntPerQuad()</code>.
%>                    @f$[N \times N \times R \times 4]@f$
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
function ret = assembleMatEdgeTetraPhiIntPhiIntFuncContNu(g, markE0T, funcCont, qOrd, refEdgePhiIntPhiIntPerQuad)
K = g.numT;
[N, ~, ~, ~] = size(refEdgePhiIntPhiIntPerQuad);
[Q, ~] = quadRule1D(qOrd); R = length(Q);
ret = { sparse(K*N, K*N), sparse(K*N, K*N) };
for n = 1 : 4
  [Q1, Q2] = gammaMapTetra(n, Q);
  funcQ0E = funcCont(g.mapRef2Phy(1, Q1, Q2), g.mapRef2Phy(2, Q1, Q2));
  for m = 1 : 2
    markAreaNuE0T = markE0T(:, n) .* g.areaE0T(:, n) .* g.nuE0T(:, n, m);
    for r = 1 : R
      ret{m} = ret{m} + kron(spdiags(markAreaNuE0T .* funcQ0E(:, r), 0, K, K), refEdgePhiIntPhiIntPerQuad(:,:,r,n));
    end % for r
  end % for m
end % for n
end % function