% Compute lumped integrals on the reference triangle, whose integrands 
% consist of all permutations of two basis functions.

%===============================================================================
%> @file ./coreintegrateRefElemPhiPhiLumped.m
%>
%> @brief Compute lumped integrals on the reference triangle, whose 
%>        integrands consist of all permutations of two basis functions.
%===============================================================================
%>
%> @brief Compute lumped integrals on the reference triangle @f$\hat{T}@f$,
%>        whose integrands consist of all permutations of two basis functions.
%>
%> It computes a matrix
%> @f$\hat{\mathsf{M}}\in\mathbb{R}^{N\times N}@f$
%> defined by
%> @f[
%> [\hat{\mathsf{M}}]_{i,j} =
%>   \sum_{l=1}^{mathrm{nEdges}} \hat{\varphi}_i (x_l) \hat{\varphi}_j(x_l) / 6.
%> @f]
%>
%> @param  N    The local number of degrees of freedom
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi2D.
%> @param  basesOnCorners  A struct containing precomputed values of the basis
%>                         functions on the corners. Must provide at
%>                         least phi2D.
%> @param  qOrd (optional) The order of the quadrature rule to be used.
%> @retval ret  The computed matrix @f$[N\times N]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2015 Florian Frank, Balthasar Reuter, Vadym Aizinger
%>
%> @author Andreas Rupp, 2018.
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
function ret = integrateRefElemPhiPhiLumped(N, basesOnQuad, basesOnCorners)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
ret = zeros(N,N); % [N x N]
for i = 1 : N
  for j = 1 : N
    for k = 1 : 3
      ret(i,j) = ret(i,j) + basesOnCorners(i,k) * basesOnCorners(j,k) / 6;
    end
  end % for
end % for
end % function
