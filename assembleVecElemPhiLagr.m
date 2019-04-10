% Assembles a vector, containing integrals of the linear Lagrange basis 
% functions.

%===============================================================================
%> @file assembleVecElemPhiLagr.m
%>
%> @brief Assembles a vector, containing integrals of the linear Lagrange  
%>        basis functions.
%===============================================================================
%>
%> @brief Assembles a vector, containing integrals of the linear Lagrange  
%>        basis functions.
%>
%> The vector @f$\mathsf{M}_L \in \mathbb{R}^{M}@f$ is component-wise by
%> @f[
%>   [\mathsf{M}_L]_{i} = \sum_{k=1}^K \int_{T_k} 
%>      \varphi_{i}^{L}\mathrm{d} \mathbf{x} \,.
%> @f]
%> where \varphi_{i}^{L} is the piecewise linear continuous function whose
%> value in the i-th global grid vertex is one and zero in all others.
%>
%> For the assembly of the local integrals we evaluate the integral by
%> the midpoint rule, resulting in the weight |T_k|/3, where |T_k| is the
%> area of element T_k.
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @retval ret        The assembled vector @f$[L \times 1]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Hennes Hajduk, Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function ret = assembleVecElemPhiLagr(g)
ret = zeros(g.numV, 1);
for k = 1 : g.numT
  for i = 1 : 3
    ret(g.V0T(k,i)) = ret(g.V0T(k,i)) + g.areaT(k) / 3;
  end % for
end % for
end % function