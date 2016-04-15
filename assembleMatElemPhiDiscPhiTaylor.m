% Assembles the global basis transformation matrix, which is required to
% transform a representation matrix between (modal) DG basis and Taylor basis.
%
%===============================================================================
%> @file assembleMatElemPhiDiscPhiTaylor.m
%>
%> @brief Assembles the global basis transformation matrix, which is required
%>        to transform a representation matrix between (modal) DG basis and 
%>        Taylor basis.
%===============================================================================
%>
%> @brief Assembles the global basis transformation matrix 
%>        @f$\mathsf{M}^\mathrm{DG,Taylor}@f$, which is required
%>        to transform a representation matrix between (modal) DG basis and 
%>        Taylor basis.
%>
%> The matrix @f$\mathsf{M}^\mathrm{DG,Taylor} \in \mathbb{R}^{KN\times KN}@f$
%> is block-diagonal and given by
%> @f[
%>  \mathsf{M}^\mathrm{DG,Taylor} = \mathrm{diag}
%>   \left( \mathsf{M}_{T_1}^\mathrm{DG,Taylor}, \ldots,
%>          \mathsf{M}_{T_K}^\mathrm{DG,Taylor} \right)
%> @f]
%> with the local transformation matrix
%> @f[
%>  \mathsf{M}_{T_k}^\mathrm{DG,Taylor} \;=\;
%>  \int_{T_k}\begin{bmatrix}
%>   \varphi_{k1}\,\phi_{k1} & \cdots & \varphi_{k1}\,\phi_{kN} ~\\
%>   \vdots                  & \ddots & \vdots \\
%>   \varphi_{kN}\,\phi_{k1} & \cdots & \varphi_{kN}\,\phi_{kN} 
%>  \end{bmatrix}\,\mathrm{d}\mathbf{x}\;,
%> @f]
%> where @f$\varphi_{ki}(\mathbf{x})@f$ denotes the ith DG basis function
%> on element @f$T_k@f$ (the basis function on the reference element is
%> given by <code>phi()</code>), and @f$\phi_{ki}@f$ is the ith Taylor
%> basis function (as provided by <code>phiTaylorRef()</code> or 
%> <code>phiTaylorPhy()</code>).
%>
%> @param  g   The lists describing the geometric and topological 
%>             properties of a triangulation (see <code>generateGridData()</code>) 
%>             @f$[1 \times 1 \text{ struct}]@f$
%> @param  N   The local number of degrees of freedom
%> @retval ret The assembled matrix @f$[KN \times KN]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function ret = assembleMatElemPhiDiscPhiTaylor(g, N, basesOnQuad)
p = (sqrt(8*N+1)-3)/2;  qOrd = max(2*p+1, 1);  [Q1,Q2,W] = quadRule2D(qOrd);
K = g.numT;
ret = sparse(K*N, K*N);
for i = 1 : N
  for j = 1 : N
    intPhiIPhiJ =  ( repmat(basesOnQuad.phi2D{qOrd}(:, i)', [K 1]) .* phiTaylorRef(g, j, Q1, Q2) ) * W';
    ret = ret + sparse(i : N : K*N, j : N : K*N, 2 * g.areaT .* intPhiIPhiJ, K*N, K*N );
  end % for
end % for
end % function