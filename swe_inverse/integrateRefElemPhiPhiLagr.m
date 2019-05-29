% Compute integrals on the reference triangle, whose integrands 
% consist of all permutations of two basis functions.

%===============================================================================
%> @file
%>
%> @brief Compute integrals on the reference triangle, whose 
%>        integrands consist of all permutations of two basis functions.
%===============================================================================
%>
%> @brief Compute integrals on the reference triangle @f$\hat{T}@f$,
%>        whose integrands consist of all permutations of two basis functions.
%>
%> It computes a matrix
%> @f$\hat{\mathsf{M}}\in\mathbb{R}^{N\times N}@f$
%> defined by
%> @f[
%> [\hat{\mathsf{M}}]_{i,j} =
%>   \int_{\hat{T}} \hat{\varphi}_i \hat{\varphi}_j \mathrm{d}\hat{\mathbf{x}} \,.
%> @f]
%>
%> @param  N    The local number of degrees of freedom
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi2D.
%> @retval ret  The computed matrix @f$[N\times N]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2015 Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function ret = integrateRefElemPhiPhiLagr(N, basesOnQuad)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
p = (sqrt(8*N+1)-3)/2;  qOrd = max(2*p, 1);  [Q1, Q2, W] = quadRule2D(qOrd);
ret = zeros(N); % [N x N]
L = [1-Q1-Q2; Q1; Q2];
for i = 1 : N
  for j = 1 : N
    ret(i, j) = sum( W .* phi(i,Q1,Q2) .* L(j,:) );
  end % for
end % for
end % function