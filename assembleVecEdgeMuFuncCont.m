% Assembles a vector containing integrals over edges of products of a edge 
% basis function with a continuous function.

%===============================================================================
%> @file
%>
%> @brief Assembles a vector containing integrals over edges of products of a 
%>        edge basis function with a continuous function.
%===============================================================================
%>
%> @brief Assembles a vector @f${\boldsymbol{K}}_{\mu,\mathrm{in}}@f$ containing integrals 
%>        over edges of products of a edge basis function with a continuous function
%>        @f$c_\mathrm{D}(t, \mathbf{x})@f$.
%> 
%> The vector @f${\boldsymbol{K}}_{\mu,\mathrm{in}} \in \mathbb{R}^{\bar{K}\bar{N}}@f$, with @f$\bar{K}@f$ being the number of edges and @f$\bar{N}@f$ being the number of edge degrees of freedom, is defined
%> component-wise by
%> @f[
%> [{\boldsymbol{K}}_{\mu,\mathrm{in}}]_{(\bar{k}-1)\bar{N}+i} = \sum_{E_{kn}\in\partial{T_{k}}\cap{\mathcal{E}}_{\text{in}}} \,
%> \int_{E_{kn}} c_\mathrm{D}(t) \, \mu_{kni}\, \text{d}s\,.
%> @f]
%> 
%> It is transformed on the unit interval as described in <code>assembleVecEdgePhiIntFuncCont()</code> which gives 
%>
%> @f[
%> \int_{E_{kn}} c_\mathrm{D}\,\mu_{kni} \, \text{d}s 
%> = {|E_{kn}|} \int_0^1 c_\mathrm{D}(t) \circ \boldsymbol{F}_{k} \circ \boldsymbol{\hat{\gamma}}_{n}(s) \, \hat{\mu}_{i}(s) \, \text{d}s \,.
%> @f]
%> 
%> and using a 1D quadrature rule provided by <code>quadRule1D()</code>. The
%> main difference is the occurence of the special edge basis function 
%> @f$\hat{\mu}@f$. We eventually obtain
%> 
%> @f[
%> \int_{E_{kn}} c_\mathrm{D}\,\mu_{kni} \, \text{d}s 
%> \approx {|E_{kn}|} \sum_{r=1}^R \omega_{r} \, c_\mathrm{D}(t) \circ \boldsymbol{F}_{k} \circ \boldsymbol{\hat{\gamma}}_{n}(\hat{q}_r) \, \hat{\mu}_{i}(\hat{q}_r) \,.
%> @f]
%> 
%> The global vector is then constructed as
%> 
%> @f[
%> {\boldsymbol{K}}_{\mu,\mathrm{in}} = 
%> \sum_{n=1}^3 {|E_{kn}|} \sum_{r=1}^R \omega_{r} \, \hat{\mu}_{i}(\hat{q}_r) \, {(\mathsf{\Delta}_n)}^\mathrm{T} \, 
%> \begin{bmatrix} \delta_{E_{1n}\in{\mathcal{E}}_{\text{in}}} & & \\ & \ddots & \\ & & \delta_{E_{Kn}\in{\mathcal{E}}_{\text{in}}} \end{bmatrix} [{\boldsymbol{C}}_\mathrm{D}]_{:,n,r} \,,
%> @f]
%> 
%> with @f$\omega_{r}@f$ denoting the integration weights, @f$\hat{q}_{r}@f$ denoting integration points, @f$\mathsf{\Delta}_n@f$ denoting the permutation matrix from element-based to edge-based view as described in <code>assembleMatEdgePhiIntMuVal()</code> and @f$\delta_{E_{kn}\in{\mathcal{E}}_{\text{in}}}@f$ denoting the Kronecker delta.
%> 
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0T    A marker indicating whether an edge should be 
%>                    recognized or not. @f$[K \times 3]@f$
%> @param  funcCont   A function handle for the continuous function
%> @param  basesOnQuad    A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least mu.
%> @param  qOrd    Order of quadrature rule to be used.
%>
%> @retval ret        The assembled vector @f$[\bar{K}\bar{N}]@f$
%>
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
%> @author Alexander Jaust, 2017
%> @author Balthasar Reuter, 2017
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
function ret = assembleVecEdgeMuFuncCont(g, markE0T, funcCont, basesOnQuad, qOrd)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
validateattributes(funcCont, {'function_handle'}, {}, mfilename, 'funcCont');
validateattributes(markE0T, {'logical'}, {'size', [g.numT 3]}, mfilename, 'markE0T');

% Determine quadrature rule
[Q, W] = quadRule1D(qOrd);
[R, Nmu] = size(basesOnQuad.mu{qOrd});

% Assemble vector
ret = zeros(g.numE, Nmu);
for n = 1 : 3
  [Q1, Q2] = gammaMap(n, Q);
  funcOnQuad = funcCont(g.mapRef2Phy(1, Q1, Q2), g.mapRef2Phy(2, Q1, Q2));
  Kkn = markE0T(:, n) .* g.areaE0T(:,n);
  ret(g.E0T(:, n), :) = ret(g.E0T(:, n), :) + (repmat(Kkn, 1, R) .* funcOnQuad) ...
                          * ( repmat(W', 1, Nmu) .* basesOnQuad.mu{qOrd} );
end % for

ret = reshape(ret', [], 1);
end % function
