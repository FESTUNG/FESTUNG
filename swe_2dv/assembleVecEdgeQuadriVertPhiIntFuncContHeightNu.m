% Assembles a vector containing integrals over vertical edges of products
% of a basis function with a continuous function, the x1-component
% of the edge normal, and divided by the smoothed mesh height.

%===============================================================================
%> @file
%>
%> @brief Assembles a vector containing integrals over vertical edges of products
%>        of a basis function with a continuous function, the @f$x^1@f$-component
%>        of the edge normal, and divided by the smoothed mesh height @f$H_s@f$.
%===============================================================================
%>
%> @brief Assembles a vector containing integrals over vertical edges of products
%>        of a basis function with a continuous function, the @f$x^1@f$-component
%>        of the edge normal, and divided by the smoothed mesh height @f$H_s@f$.
%>
%>
%> The vector @f$\mathbf{J}_\mathrm{D} \in \mathbb{R}^{KN}@f$ is defined
%> component-wise by
%> @f[
%> [\mathbf{J}_\mathrm{D}]_{(k-1)N+i} =
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}^v_D}
%>  \nu_{kn}^1 \int_{E_{kn}} \frac{1}{H_s} \varphi_{ki} c_\mathrm{D}(t) \mathrm{d}s\,,
%> @f]
%> with @f$\nu_{kn}^1@f$ the @f$x^1@f$-component of the edge normal.
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
%> @param  funcCont   A function handle for the continuous function @f$c_\mathrm{D}@f$.
%> @param heightV0T1D The smoothed mesh height in each node of the 1D grid
%>                    @f$[\overline{K} \times 2]@f$
%> @param  N          The number of local degrees of freedom @f$[\text{scalar}]@f$
%> @param  qOrd       The order of the quadrature rule to be used. 
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                    functions on quadrature points. Must provide at
%>                    least phi1D.
%>
%> @retval ret        The assembled vectors @f$[KN \times 1]@f$
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
function ret = assembleVecEdgeQuadriVertPhiIntFuncContHeightNu(g2D, g1D, markE0T, funcCont, heightV0T1D, N, qOrd, basesOnQuad)
K = g2D.numT;
[Q, W] = quadRule1D(qOrd);
ret = zeros(K, N);
for n = 3 : 4
  [Q1, Q2] = gammaMapQuadri(n, Q);
  funcQ0E = funcCont(g2D.mapRef2Phy(1, Q1, Q2), g2D.mapRef2Phy(2, Q1, Q2));
  areaNuE0T = markE0T(:,n) .* g2D.areaE0T(:,n) .* g2D.nuE0T(:,n,1);
  areaNuE0THeight = areaNuE0T ./ (g1D.markT2DT * heightV0T1D(:,5-n));
  for i = 1 : N
    ret(:, i) = ret(:, i) + areaNuE0THeight .* ( funcQ0E * (W.' .* basesOnQuad.phi1D{qOrd}(:, i, n)) );
  end  % for i
end  % for n
ret = reshape(ret.', K*N, 1);
end  % function