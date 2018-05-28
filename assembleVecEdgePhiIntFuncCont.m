% Assembles a vector containing integrals over edges of products of a basis 
% function with a continuous function.

%===============================================================================
%> @file
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
%>  \int_{E_{kn}} \varphi_{ki} c_\mathrm{D}(t) \mathrm{d}s\,.
%> @f]
%>
%> For the implementation, the integrals are backtransformed to the
%> reference element @f$\hat{T}@f$ and further to the unit interval @f$[0,1]@f$
%> using the mapping @f$\hat{\mathbf{\gamma}}_n(s)@f$ as provided by 
%> <code>gammaMap()</code>.
%> This gives the component-wise formulation
%> @f[
%>  [\mathbf{K}_\mathrm{D}]_{(k-1)N+i} =
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_D} |E_{kn}|
%>  \int_0^1 \hat{\varphi}_{i} \circ \hat{\mathbf{\gamma}}_n(s)
%>     c_\mathrm{D}(t, \mathbf{F}_k \circ \hat{\mathbf{\gamma}}_n(s))
%>     \mathrm{d}s \,.
%> @f]
%> This integral is then approximated using a 1D quadrature rule provided by
%> <code>quadRule1D()</code> 
%> @f[
%>  [\mathbf{K}_\mathrm{D}]_{(k-1)N+i} \approx
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_D} |E_{kn}|
%>  \sum_{r=1}^R \omega_r \hat{\varphi}_{i} \circ \hat{\mathbf{\gamma}}_n(q_r)
%>     c_\mathrm{D}(t, \mathbf{F}_k \circ \hat{\mathbf{\gamma}}_n(q_r)) \,,
%> @f]
%> allowing to vectorize over all elements.
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0T    <code>logical</code> arrays that mark each elements
%>                    edges on which the vector entries should be
%>                    assembled @f$[K \times n_\mathrm{edges}]@f$
%> @param  funcCont   A function handle for the continuous function
%> @param  N          The number of local degrees of freedom @f$[\text{scalar}]@f$
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi1D.
%> @param qOrd        (optional) Order of quadrature rule to be used.
%                     Defaults to @f$2p+1@f$.
%> @param coefE0T     (optional) Coefficient vector that is applied to each
%>                    block. Defaults to <code>g.areaE0T</code>
%>                    @f$[K \times n_\mathrm{edges}]@f$
%> @retval ret        The assembled vector @f$[KN]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Florian Frank, Balthasar Reuter, Vadym Aizinger
%>
%> @author Balthasar Reuter, 2017.
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
function ret = assembleVecEdgePhiIntFuncCont(g, markE0T, funcCont, N, basesOnQuad, qOrd, coefE0T)
% Extract dimensions
K = g.numT;  nEdges = size(g.E0T, 2);

% Determine quadrature rule and mapping to physical element
if nargin < 6
  switch nEdges
    case 3
      p = (sqrt(8*N+1)-3)/2;
    case 4
      p = sqrt(N)-1;
    otherwise
      error('Unknown number of edges')
  end % switch
  qOrd = 2*p+1;  
end % if
[Q, W] = quadRule1D(qOrd);

% Determine default coefficient
if nargin < 7
  coefE0T = g.areaE0T;
end % if

% Check function arguments that are directly used
validateattributes(funcCont, {'function_handle'}, {}, mfilename, 'funcCont');
validateattributes(markE0T, {'logical'}, {'size', [K nEdges]}, mfilename, 'markE0T');
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
validateattributes(coefE0T, {'numeric'}, {'size', [K nEdges]}, mfilename, 'coefE0T');

% Determine edge mapping
switch nEdges
  case 3
    gamma = @gammaMap;
  case 4
    gamma = @gammaMapQuadri;
  otherwise
    error('Unknown number of edges')
end % switch

% Assemble vector
ret = zeros(K, N);
for n = 1 : nEdges
  [Q1, Q2] = gamma(n, Q);
  funcQ0E = funcCont(g.mapRef2Phy(1, Q1, Q2), g.mapRef2Phy(2, Q1, Q2));
  markCoefE0T = markE0T(:, n) .* coefE0T(:, n);
  for i = 1 : N
    ret(:,i) = ret(:,i) + markCoefE0T .* ( funcQ0E * (W' .* basesOnQuad.phi1D{qOrd}(:, i, n)) );
  end % for
end % for
ret = reshape(ret.', K*N, 1);
end % function
