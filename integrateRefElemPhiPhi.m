% Compute integrals on the reference triangle, whose integrands 
% consist of all permutations of two basis functions.
%
%===============================================================================
%> @file integrateRefElemPhiPhi.m
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
%> @param  N    The local number of degrees of freedom, either as a scalar
%>              for both basis functions, or as a vector with two
%>              entries, specifying the number of degrees of freedom for
%>              each basis function.
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi2D.
%> @retval ret  The computed matrix @f$[N\times N]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Florian Frank, Balthasar Reuter, Vadym Aizinger
%>											Modified by Hennes Hajduk 06/30/2016
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
function ret = integrateRefElemPhiPhi(N, basesOnQuad)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
p = (sqrt(8*max(N)+1)-3)/2;  qOrd = max(2*p, 1);  [~,~,W] = quadRule2D(qOrd);

if length(N) == 1
  N = N * ones(2,1);
else
  validateattributes(N, {'numeric'}, {'numel', 2}, mfilename, 'N')
end % if

ret = zeros(N(1), N(2));
for i = 1 : N(1)
  for j = 1 : N(2)
    ret(i, j) = sum( W' .* basesOnQuad.phi2D{qOrd}(:, i) .* basesOnQuad.phi2D{qOrd}(:, j) );
  end % for
end % for
end % function
