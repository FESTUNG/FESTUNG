% Compute integrals on the reference interval, whose integrands consist of all
% permutations of a basis function with an edge basis function.
%
%===============================================================================
%> @file integrateRefEdgePhiIntMu.m
%>
%> @brief NEW Compute integrals on the reference interval, whose integrands 
%> consist of all permutations of a basis function with an edge basis function.
%===============================================================================
%>
%> @brief Compute integrals on the reference interval @f$[0,1]@f$, whose 
%> integrands consist of all permutations of a basis function with an edge 
%> basis function.
%>
%> It computes a multidimensional array
%> @f$\mathsf{\hat{R}}_{\mu}\in\mathbb{R}^{N\times \bar{N}\times 3 \times 2}@f$
%> defined by
%> @f[
%> [\mathsf{\hat{R}}_{\mu}]_{i,j,n,l} := \int_{0}^{1} \hat{\varphi}_{i} \circ \boldsymbol{\hat{\gamma}}_{n}(s) \, \hat{\mu}_{j} \circ \hat{\beta}_{kn}(s) \, \text{d}s\, ,
%> @f]
%> where the mapping @f$\hat{\mathbf{\gamma}}_n@f$ is given in 
%> <code>gammaMap()</code> and @f$\hat{\beta}_{kn} @f$the mapping as described in <code>TODO</code>.
%> @param  N    The local number of degrees of freedom @f$[2 \text{vector}]@f$
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi1D, mu and thetaMu.  @f$[\text{struct}]@f$
%> @retval ret  The computed array @f$[N\times \bar{N}\times 3 \times 2]@f$
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
function ret = integrateRefEdgePhiIntMu(N, basesOnQuad, qOrd)
validateattributes(N, {'numeric'}, {'numel', 2}, mfilename, 'N')
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
if nargin < 3, p = (sqrt(8*N(1)+1)-3)/2;  qOrd = 2*p+1;  end
[~, W] = quadRule1D(qOrd);

ret = zeros(N(1), N(2), 3, 2); % [N(1) x N(2) x 3 x 2]
for n = 1 : 3
    for i = 1 : N(1)
        for j = 1 : N(2)
            ret(i, j, n, 1) =  W * (basesOnQuad.phi1D{qOrd}(:, i, n) .* basesOnQuad.mu{qOrd}(:, j) );
            ret(i, j, n, 2) =  W * (basesOnQuad.phi1D{qOrd}(:, i, n) .* basesOnQuad.thetaMu{qOrd}(:, j, 2) );
        end % for
    end % for
end % for
end % function
