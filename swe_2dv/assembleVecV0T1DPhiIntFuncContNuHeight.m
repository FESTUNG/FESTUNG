% Assembles a vector containing products of a 1D basis function with a 
% continuous depth-integrated function, evaluated at the endpoints of the 1D 
% element, the x1-component of the "normal", and divided by the smoothed mesh height.

%===============================================================================
%> @file
%>
%> @brief Assembles a vector containing products of a 1D basis function with a 
%>        continuous depth-integrated function, evaluated at the endpoints of the
%>        1D element, the x1-component of the "normal", and divided by the 
%>        smoothed mesh height.
%===============================================================================
%>
%> @brief Assembles a vector containing products of a 1D basis function with a 
%>        continuous depth-integrated function, evaluated at the endpoints of the
%>        1D element, the x1-component of the "normal", and divided by the 
%>        smoothed mesh height.
%>
%>
%> The vector @f$\overline{\vec{J}}_\mathrm{D} \in \mathbb{R}^{\overline{KN}}@f$ is defined
%> component-wise by
%> @f[
%> \left[\overline{\vec{J}}_\mathrm{D}{D}\right]_{(\overline{k}-1)\overline{N}+i} :=
%> \sum_{E_{kn}\in\Pi^{-1}\partial \overline{T}_{\overline{k}} \cap \mathcal{E}_\mathrm{D}}
%> \nu_{kn}^1 \int_{E_{kn}} \frac{1}{H_s} \,\phi_{\overline{k}i}  \, 
%> u_\mathrm{D} \, \mathrm{d}\sigma \,,
%> @f]
%> with @f$\nu_{\overline{k}\overline{n}}^1@f$ the @f$x^1@f$-component of the edge
%> normal (i.e., @f$\pm 1@f$) and @f$u_\mathrm{D}@f$ the @f$(x^1,x^2)@f$-dependent
%> continuous function and @f$H_s@f$ the smoothed mesh height.
%>
%> With @f$\overline{u}_\mathrm{D} := \int_{\zeta_b}^\xi u_\mathrm{D} \mathrm{d}x^2@f$
%> the depth-integrated function, this simplifies to an evaluation of all 
%> quantities in the endpoints of the 1D element.
%>
%> @param  g2D        The lists describing the geometric and topological 
%>                    properties of a quadrilateral triangulation (see 
%>                    <code>domainRectTrap()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  g1D        The lists describing the geometric and topological 
%>                    properties of a matching 1D triangulation (see 
%>                    <code>generateGridData1D()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0T    <code>logical</code> array that marks each elements
%>                    edges on which the vector entries should be
%>                    assembled @f$[K \times 4]@f$
%> @param  funcCont   A function handle for the continuous function @f$u_\mathrm{D}@f$.
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
function ret = assembleVecV0T1DPhiIntFuncContNuHeight(g2D, g1D, markE0T, funcCont, heightV0T1D, barN, qOrd, basesOnQuad)
barK = g1D.numT;
[Q, W] = quadRule1D(qOrd);
ret = zeros(barK, barN);
for n = 3 : 4
  [Q1, Q2] = gammaMapQuadri(n, Q);
  funcQ0E = funcCont(g2D.mapRef2Phy(1, Q1, Q2), g2D.mapRef2Phy(2, Q1, Q2));
  areaNuE0T = markE0T(:,n) .* g2D.areaE0T(:,n) .* g2D.nuE0T(:,n,1);
  integrand = g1D.markT2DT.' * (areaNuE0T .* (funcQ0E * W.'));
  ret = ret + (integrand ./ heightV0T1D(:, 5-n)) * basesOnQuad.phi0D{qOrd}(:, 5-n).';
end  % for n
ret = reshape(ret.', barK*barN, 1);
end  % function