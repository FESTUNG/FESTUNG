% Assembles a vector containing integrals over edges of products of a basis 
% function with a continuous function.

%===============================================================================
%> @file assembleVecEdgeTetraPhiIntFuncCont.m
%>
%> @brief Assembles a vector containing integrals over edges of products of a 
%>        basis function with a continuous function.
%===============================================================================
%>
%> @brief Assembles a vector @f$\mathbf{K}_\mathrm{D}@f$ containing integrals 
%>        over edges of products of a basis function with a continuous function
%>        @f$c_\mathrm{D}(t, \mathbf{x})@f$.
%>
%> The vector @f$\mathbf{K}_\mathrm{D} \in \mathbb{R}^{KN}@f$ is defined
%> component-wise by
%> @f[
%> [\mathbf{K}_\mathrm{D}]_{(k-1)N+i} =
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_D}
%>  \frac{1}{|E_{kn}|} \int_{E_{kn}} \varphi_{ki} c_\mathrm{D}(t) \mathrm{d}s\,.
%> @f]
%> For the implementation, the integrals are backtransformed to the
%> reference square @f$\hat{T} = [0,1]^2@f$ using a mapping
%> @f$\mathbf{F}_k:\hat{T}\ni\hat{\mathbf{x}}\mapsto\mathbf{x}\in T_k@f$
%> This allows to reformulate 
%> @f[
%>   \int_{E_{kn}} \varphi_{ki} c_\mathrm{D}(t) \mathrm{d}s =
%>   \frac{|E_{kn}|}{|\hat{E}_n|} \int_{\hat{E}_n} \hat{\varphi}_i 
%>   c_\mathrm{D}(t,\mathbf{F}_k(\hat{\mathbf{x}}))\mathrm{d}\hat{\mathbf{x}}\,.
%> @f]
%> Further transformation to the unit interval @f$[0,1]@f$ using the mapping 
%> @f$\hat{\mathbf{\gamma}}_n(s)@f$ as provided by <code>gammaMapTetra()</code>
%> gives the component-wise formulation
%> @f[
%>  [\mathbf{K}_\mathrm{D}]_{(k-1)N+i} =
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_D}
%>  \int_0^1 \hat{\varphi}_{i} \circ \hat{\mathbf{\gamma}}_n(s)
%>     c_\mathrm{D}(t, \mathbf{F}_k \circ \hat{\mathbf{\gamma}}_n(s))
%>     \mathrm{d}s \,.
%> @f]
%> This integral is then approximated using a 1D quadrature rule provided by
%> <code>quadRule1D()</code> 
%> @f[
%>  [\mathbf{K}_\mathrm{D}]_{(k-1)N+i} \approx
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_D}
%>  \sum_{r=1}^R \omega_r \hat{\varphi}_{i} \circ \hat{\mathbf{\gamma}}_n(q_r)
%>     c_\mathrm{D}(t, \mathbf{F}_k \circ \hat{\mathbf{\gamma}}_n(q_r)) \,,
%> @f]
%> allowing to vectorize over all triangles.
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>domainRectTrap()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0T    <code>logical</code> arrays that mark each elements
%>                    edges on which the vector entries should be
%>                    assembled @f$[K \times 3]@f$
%> @param  funcCont   A function handle for the continuous function
%> @param  N          The number of local degrees of freedom @f$[\text{scalar}]@f$
%> @param  qOrd       The order of the quadrature rule to be used.
%> @param coefE0T     (optional) Coefficient vector that is applied to each
%>                    block. Defaults to <code>g.areaE0T</code>
%>                    @f$[K \times 4]@f$
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi1D.
%> @retval ret        The assembled vector @f$[KN]@f$
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
function ret = assembleVecEdgeTetraPhiIntFuncCont(g, markE0T, funcCont, N, qOrd, basesOnQuad, coefE0T)
if nargin < 7
  coefE0T = g.areaE0T;
end % if
K = g.numT;
[Q, W] = quadRule1D(qOrd);
ret = zeros(K, N);
for n = 1 : 4
  [Q1, Q2] = gammaMapTetra(n, Q);
  funcQ0E = funcCont(g.mapRef2Phy(1, Q1, Q2), g.mapRef2Phy(2, Q1, Q2));
  markCoefE0T = markE0T(:, n) .* coefE0T(:, n);
  for i = 1 : N
    ret(:, i) = ret(:, i) + markCoefE0T .* ( funcQ0E * ( W.' .* basesOnQuad.phi1D{qOrd}(:, i, n) ) );
  end  % for i
end  % for n
ret = reshape(ret.', K*N, 1);
end  % function