% Assembles a matrix, containing integrals of products of two Taylor 
% basis functions.

%===============================================================================
%> @file assembleMatElemPhiTaylorPhiTaylor.m
%>
%> @brief Assembles a matrix, containing integrals of products of two Taylor 
%>        basis functions. This corresponds to a mass matrix.
%===============================================================================
%>
%> @brief Assembles a mass matrix @f$\mathsf{M}^\mathrm{Taylor}@f$
%>        containing integrals of products of two Taylor basis functions.
%>
%> The matrix @f$\mathsf{M}^\mathrm{Taylor} \in \mathbb{R}^{KN\times KN}@f$ 
%> is block diagonal and defined component-wise by
%> @f[
%>   [\mathsf{M}^\mathrm{Taylor}]_{(k-1)N+i,(k-1)N+j} = \int_{T_k} 
%>      \phi_{ki} \phi_{kj} \mathrm{d} \mathbf{x} \,.
%> @f]
%> All other entries are zero.
%> The basis functions are evaluated using <code>phiTaylorRef()</code>.
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  N          The number of local degrees of freedom.
%> @retval ret        The assembled matrix @f$[KN \times KN]@f$
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
function ret = assembleMatElemPhiTaylorPhiTaylor(g, N)
% Extract dimensions and determine quadrature rule
p = (sqrt(8*N+1)-3)/2;  qOrd = max(2*p, 1); 
[Q1, Q2, W] = quadRule2D(qOrd);
K = g.numT; 

% Assemble matrix
ret = sparse(K*N, K*N);
for i = 1 : N
  for j = 1 : N
    intPhiIPhiJ =  ( phiTaylorRef(g, i, Q1, Q2) .* phiTaylorRef(g, j, Q1, Q2) ) * W';
    ret = ret + sparse(i : N : K*N, j : N : K*N, 2 * g.areaT .* intPhiIPhiJ, K*N, K*N );
  end % for
end % for
end % function
