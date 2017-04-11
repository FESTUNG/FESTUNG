% Assembles two vectors containing integrals over edges of products of a basis 
% function with a continuous function and one of the components of the edge 
% normal.

%===============================================================================
%> @file assembleVecEdgeTetraPhiIntFuncContNu.m
%>
%> @brief Assembles two vectors containing integrals over edges of products of a 
%>        basis function with a continuous function and one of the components of
%>        the edge normal.
%===============================================================================
%>
%> @brief Assembles two vectors @f$\mathbf{J}_\mathrm{D}^m, m\in\{1,2\}@f$ 
%>        containing integrals over edges of products of a basis function with a
%>        continuous function @f$c_\mathrm{D}(t, \mathbf{x})@f$ and one of the 
%>        components of the edge normal.
%>
%> The vectors @f$\mathbf{J}_\mathrm{D}^m \in \mathbb{R}^{KN}@f$ are defined
%> component-wise by
%> @f[
%> [\mathbf{J}_\mathrm{D}^m]_{(k-1)N+i} =
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_D}
%>  \nu_{kn}^m \int_{E_{kn}} \varphi_{ki} c_\mathrm{D}(t) \mathrm{d}s\,,
%> @f]
%> with @f$\nu_{kn}@f$ the @f$m@f$-th component (@f$m\in\{1,2\}@f$) of the edge
%> normal.
%>
%> For the implementation, the integrals are backtransformed to the
%> reference square @f$\hat{T} = [0,1]^2@f$ using a mapping
%> @f$\mathbf{F}_k:\hat{T}\ni\hat{\mathbf{x}}\mapsto\mathbf{x}\in T_k@f$.
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
%>  [\mathbf{J}_\mathrm{D}^m]_{(k-1)N+i} =
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_D} \nu_{kn}^m |E_{kn}|
%>  \int_0^1 \hat{\varphi}_{i} \circ \hat{\mathbf{\gamma}}_n(s)
%>     c_\mathrm{D}(t, \mathbf{F}_k \circ \hat{\mathbf{\gamma}}_n(s))
%>     \mathrm{d}s \,.
%> @f]
%> This integral is then approximated using a 1D quadrature rule provided by
%> <code>quadRule1D()</code> 
%> @f[
%>  [\mathbf{J}_\mathrm{D}^m]_{(k-1)N+i} \approx
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_D} \nu_{kn}^m |E_{kn}|
%>  \sum_{r=1}^R \omega_r \hat{\varphi}_{i} \circ \hat{\mathbf{\gamma}}_n(q_r)
%>     c_\mathrm{D}(t, \mathbf{F}_k \circ \hat{\mathbf{\gamma}}_n(q_r)) \,,
%> @f]
%> allowing to vectorize over all triangles.
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0T    <code>logical</code> arrays that mark each triangles
%>                    edges on which the vector entries should be
%>                    assembled @f$[K \times 4]@f$
%> @param  funcCont   A function handle for the continuous function
%> @param  N          The number of local degrees of freedom @f$[\text{scalar}]@f$
%> @param  qOrd       The order of the quadrature rule to be used.
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi1D.
%> @retval ret        The assembled vectors @f$2\times1 \text{ cell}@f$
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
function ret = assembleVecEdgeTetraPhiIntFuncContNu(g, markE0T, funcCont, N, qOrd, basesOnQuad)
if iscell(funcCont)
  ret = assembleVecEdgeTetraPhiIntFuncContNuVector(g, markE0T, funcCont, N, qOrd, basesOnQuad);
else
  ret = assembleVecEdgeTetraPhiIntFuncContNuScalar(g, markE0T, funcCont, N, qOrd, basesOnQuad);
end % if
end % function

function ret = assembleVecEdgeTetraPhiIntFuncContNuVector(g, markE0T, funcCont, N, qOrd, basesOnQuad)
K = g.numT;
[Q, W] = quadRule1D(qOrd);
ret = { zeros(K, N), zeros(K, N) };
for n = 1 : 4
  [Q1, Q2] = gammaMapTetra(n, Q);
  for m = 1 : 2
    funcQ0E = funcCont{1}(g.mapRef2Phy(1, Q1, Q2), g.mapRef2Phy(2, Q1, Q2));
    intQ0E = funcQ0E * (repmat(W(:), 1, N) .* basesOnQuad.phi1D{qOrd}(:, :, n));
    ret{m} = ret{m} + bsxfun(@times, markE0T(:, n) .* g.areaE0T(:, n) .* g.nuE0T(:, n, m), intQ0E);
  end  % for m
end  % for n
ret = cellfun(@(c) reshape(c.', K*N, 1), ret, 'UniformOutput', false);
end % function

function ret = assembleVecEdgeTetraPhiIntFuncContNuScalar(g, markE0T, funcCont, N, qOrd, basesOnQuad)
K = g.numT;
[Q, W] = quadRule1D(qOrd);
ret = { zeros(K, N), zeros(K, N) };
for n = 1 : 4
  [Q1, Q2] = gammaMapTetra(n, Q);
  funcQ0E = funcCont(g.mapRef2Phy(1, Q1, Q2), g.mapRef2Phy(2, Q1, Q2));
  intQ0E = funcQ0E * (repmat(W(:), 1, N) .* basesOnQuad.phi1D{qOrd}(:, :, n));
  for m = 1 : 2
    ret{m} = ret{m} + bsxfun(@times, markE0T(:, n) .* g.areaE0T(:, n) .* g.nuE0T(:, n, m), intQ0E);
  end % for m
end  % for n
ret = cellfun(@(c) reshape(c.', K*N, 1), ret, 'UniformOutput', false);
end  % function