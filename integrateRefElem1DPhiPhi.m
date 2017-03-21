% Compute integrals on the reference interval, whose integrands 
% consist of all permutations of two basis functions.

%===============================================================================
%> @file integrateRefElem1DPhiPhi.m
%>
%> @brief Compute integrals on the reference interval, whose 
%>        integrands consist of all permutations of two basis functions.
%===============================================================================
%>
%> @brief Compute integrals on the reference interval, whose 
%>        whose integrands consist of all permutations of two basis functions.
%>
%> It computes a matrix
%> @f$\hat{\overline{\mathsf{M}}}\in\mathbb{R}^{N\times N}@f$
%> defined by
%> @f[
%> [\hat{\overline{\mathsf{M}}}]_{i,j} =
%>   \int_0^1 \hat{\Phi}_i \hat{\Phi}_j \mathrm{d}\hat{\mathbf{x}} \,.
%> @f]
%>
%> @param  N    The local number of degrees of freedom
%> @param  qOrd The order of the quadrature rule to be used
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi1D.
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
function ret = integrateRefElem1DPhiPhi(N, qOrd, basesOnQuad)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
[~, W] = quadRule1D(qOrd);
ret = zeros(N); % [N x N]
for i = 1 : N
  for j = 1 : N
    ret(i, j) = sum( W' .* basesOnQuad.phi1D{qOrd}(:, i) .* basesOnQuad.phi1D{qOrd}(:, j) );
  end % for
end % for
end % function
