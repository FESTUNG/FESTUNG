% Assembles a vector containing products of a 1D basis function with a 
% discrete depth-integrated function, evaluated at the endpoints of the 1D 
% element, the x1-component of the "normal", and divided by the smoothed mesh height.

%===============================================================================
%> @file
%>
%> @brief Assembles a vector containing products of a 1D basis function with a 
%>        discrete depth-integrated function, evaluated at the endpoints of the
%>        1D element, the x1-component of the "normal", and divided by the 
%>        smoothed mesh height.
%===============================================================================
%>
%> @brief Assembles a vector containing products of a 1D basis function with a 
%>        discrete depth-integrated function, evaluated at the endpoints of the
%>        1D element, the x1-component of the "normal", and divided by the 
%>        smoothed mesh height.
%>
%>
%> The vector @f$\overline{\vec{J}}_\mathrm{D} \in \mathbb{R}^{\overline{KN}}@f$ is defined
%> component-wise by
%> @f[
%> \left[\overline{\vec{J}}_\mathrm{D}{D}\right]_{(\overline{k}-1)\overline{N}+i} :=
%> \sum_{a^1_{\overline{k}\overline{n}}\in\partial \overline{T}_{\overline{k}}
%>    \cap \mathcal{V}_\mathrm{D}}
%> \frac{\nu_{\overline{k}\overline{n}}^1}{H_s\left(a^1_{\overline{k}\overline{n}}\right)}\,
%> \phi_{\overline{k}i} \left(a^1_{\overline{k}\overline{n}}\right) \, 
%> \left(\sum_{j=1}^{\overline{N}} \overline{U}_{\overline{k}j}\left(a^1_{\overline{k}\overline{n}}\right)
%>  \,\phi_{\overline{k}j}\left(a^1_{\overline{k}\overline{n}}\right)\right) \,,
%> @f]
%> with @f$\nu_{\overline{k}\overline{n}}^1@f$ the @f$x^1@f$-component of the edge
%> normal (i.e., @f$\pm 1@f$) and 
%> @f$u_\Delta = \sum_{j=1}^{\overline{N}} \overline{U}_{\overline{k}j}
%>   \left(a^1_{\overline{k}\overline{n}}\right) \,\phi_{\overline{k}j}
%>   \left(a^1_{\overline{k}\overline{n}}\right)@f$ 
%  discrete function and @f$H_s@f$ the smoothed mesh height.
%>
%> @param  g1D        The lists describing the geometric and topological 
%>                    properties of a 1D triangulation (see 
%>                    <code>generateGridData1D()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markV0T    <code>logical</code> array that marks each 1D elements
%>                    nodes on which the vector entries should be
%>                    assembled @f$[\overline{K} \times 2]@f$
%> @param dataDisc    Representation matrices of the discrete function 
%>                    @f$\overline{u}_\Delta(x^1)@f$
%>                    @f$[2 \times 1 \text{ cell}]@f$
%> @param heightV0T1D The smoothed mesh height in each node of the 1D grid
%>                    @f$[\overline{K} \times 2]@f$
%> @param  barN       The number of local degrees of freedom @f$[\text{scalar}]@f$
%> @param  qOrd       The order of the quadrature rule to be used. 
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                    functions on quadrature points. Must provide at
%>                    least phi0D.
%>
%> @retval ret        The assembled vector @f$[\overline{KN} \times 1]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2018 Balthasar Reuter, Florian Frank, Vadym Aizinger
%>
%> @author Balthasar Reuter, 2018.
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
function ret = assembleVecV0T1DPhiIntFuncDiscIntNuHeight(g1D, markV0T, dataDisc, heightV0T1D, barN, qOrd, basesOnQuad)
barK = g1D.numT;
ret = zeros(barK, barN);
for n = 1 : 2
  funcDiscV0T = (dataDisc{1} + (n-1) * dataDisc{2}) * basesOnQuad.phi0D{qOrd}(:, n);
  ret = ret + (markV0T(:, n) .* g1D.nuV0T(:, n) .* funcDiscV0T ./ heightV0T1D(:, n)) * basesOnQuad.phi0D{qOrd}(:, n).';
end  % for n
ret = reshape(ret.', barK*barN, 1);
end  % function