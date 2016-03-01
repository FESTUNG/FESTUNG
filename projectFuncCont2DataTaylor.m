% Compute the Taylor basis representation of an algebraic function.
%
%===============================================================================
%> @file projectFuncCont2DataTaylor.m
%>
%> @brief Compute the Taylor basis representation of an algebraic function.
%===============================================================================
%>
%> @brief Compute the DG/modal basis representation of an algebraic function by
%>        performing the @f$L^2@f$-projection.
%>
%> The Taylor basis representation of an algebraic function @f$d(t)@f$,
%> is given as
%> @f[
%>   d_h(t, \mathbf{x}) = \sum_{j=1}^N D_{kj}(t) \phi_{kj}(\mathbf{x}) \,,
%> @f]
%> such that @f$d_h(t) \in \mathbb{P}_d (\mathcal{T}_h)@f$ is an adequate
%> approximation of @f$d(t)@f$.
%>
%> To produce @f$d_h(t)@f$ we use the @f$L^2@f$-projection defined locally for
%> @f$T_k \in \mathcal{T}_h@f$ by
%> @f[
%>  \forall w_h \in \mathbb{P}_d(T_k), \quad
%>  \int_{T_k} w_h d_h(t) = \int_{T_k} w_h d(t) \,.
%> @f]
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  funcCont   A function handle for the continuous function
%> @param  ord        The order of the quadrature rule provided by 
%>                    <code>quadRule2D()</code>
%> @param  globM      Global mass matrix @f$\mathsf{M}^\mathrm{Taylor}@f$ 
%>                    as provided by 
%>                    <code>assembleMatElemPhiTaylorPhiTaylor()</code>.
%>                    @f$[N \times N]@f$
%> @retval The representation matrix of the Taylor basis representation.
%>         @f$[K \times N]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2015 Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function dataDisc = projectFuncCont2DataTaylor(g, funcCont, ord, globM)
assert(isa(funcCont, 'function_handle'), 'funcCont must be a function_handle')
ord = max(ord,1);  [Q1, Q2, W] = quadRule2D(ord);
K = g.numT; N = size(globM, 1) / K;
F1 = @(X1, X2) g.B(:,1,1)*X1 + g.B(:,1,2)*X2 + g.coordV0T(:,1,1)*ones(size(X1));
F2 = @(X1, X2) g.B(:,2,1)*X1 + g.B(:,2,2)*X2 + g.coordV0T(:,1,2)*ones(size(X1));
rhs = zeros(K, N);
for i = 1 : N
  funcOnQuad = funcCont(F1(Q1, Q2), F2(Q1, Q2));
  phiOnQuad = phiTaylorPhy(g, i, F1(Q1, Q2), F2(Q1, Q2), ord);
  rhs(:, i) = g.areaT .* (( funcOnQuad .* phiOnQuad ) * (2 * W)');
end
dataDisc = reshape(globM \ reshape(rhs', [K*N 1]), [N K])';
end % function
