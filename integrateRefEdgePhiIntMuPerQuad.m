% Compute integral contributions on the reference interval, whose integrands 
% consist of all permutations of a basis function with an edge basis function 
% at every integration point.

%===============================================================================
%> @file
%>
%> @brief Compute integral contributions on the reference interval, whose 
%> integrands consist of all permutations of a basis function with an edge basis
%> function at every integration point.
%===============================================================================
%>
%> @brief Compute integral contributions on the reference interval @f$[0,1]@f$, whose 
%> integrands consist of all permutations of a basis function with an edge basis
%> function at every integration point.
%>
%> It computes a multidimensional array
%> @f$\mathsf{\hat{S}}\in \mathbb{R}^{ N \times \bar{N} \times 3 \times R \times 2}@f$
%> defined by
%> @f[
%> [\mathsf{\hat{S}}]_{i,j,n,r,l} = \omega_{r} \,  \hat{\varphi}_{i} \circ \boldsymbol{\hat{\gamma}}_{n}(\hat{q}_r) \, \hat{\mu}_{j}(\hat{q}_r \circ \hat{\beta}_{kn}(\hat{q}_r) 
%> @f]
%> where the mapping @f$\hat{\mathbf{\gamma}}_n@f$ is given in 
%> <code>gammaMap()</code>, @f$\omega_r@f$ is the integration weight associated with integration points @f$\hat{q}_r@f$and @f$\hat{\beta}_{kn} @f$the mapping as described in computeBasesOnQuadEdge(). For efficient access of the array entries it is actually stored as an @f$2 \times 1 \text{ cell}@f$ with each cell storing a @f$N \times \bar{N} \times 3 \times R@f$ array.
%> 
%> @param  N    The local number of degrees of freedom @f$[2 \text{ vector}]@f$
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi1D, mu and thetaMu. @f$[\text{struct}]@f$
%> @param  qOrd 	(optional) The order of the quadrature rule to be used. @f$[\text{scalar}]@f$
%>
%> @retval ret  The computed array @f$2 \times 1 \text{ cell}@f$ 
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
function ret = integrateRefEdgePhiIntMuPerQuad(N, basesOnQuad, qOrd)
validateattributes(N, {'numeric'}, {'numel', 2}, mfilename, 'N')
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
if nargin < 3, p = (sqrt(8*N+1)-3)/2;  qOrd = 2*p+1; end
[~, W] = quadRule1D(qOrd); R = length(W);
ret = { zeros(N(1), N(2), 3, R); zeros(N(1), N(2), 3, R) };
for n = 1 : 3 
  for i = 1 : N(1)
    for j = 1 : N(2)
      ret{1}(i, j, n, :) = W(:) .* basesOnQuad.phi1D{qOrd}(:, i, n) ...
                              .* basesOnQuad.mu{qOrd}(:, j);
      ret{2}(i, j, n, :) = W(:) .* basesOnQuad.phi1D{qOrd}(:, i, n) ...
                              .* basesOnQuad.thetaMu{qOrd}(:, j, 2);
    end % for j
  end % for i
end % for n
end % function
