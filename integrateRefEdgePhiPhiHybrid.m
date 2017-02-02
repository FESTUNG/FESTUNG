% Compute integrals over the edges of the reference triangle, whose integrands 
% consist of all permutations of two basis functions.

%===============================================================================
%> @file integrateRefEdgePhiIntPhiInt.m
%>
%> @brief Compute integrals over the edges of the reference triangle, whose 
%>        integrands consist of all permutations of two basis functions.
%===============================================================================
%>
%> @brief Compute integrals over the edges of the reference triangle 
%>        @f$\hat{T}@f$, whose integrands consist of all permutations of two
%>        basis functions.
%>
%> It computes a multidimensional array
%> @f$\hat{\mathsf{{S}}}^\mathrm{diag}\in\mathbb{R}^{N\times N\times3}@f$
%> defined by
%> @f[
%> [\hat{\mathsf{{S}}}^\mathrm{diag}]_{i,j,n} =
%>   \int_0^1 \hat{\varphi}_i \circ \hat{\mathbf{\gamma}}_n(s) 
%>   \hat{\varphi}_j\circ \hat{\mathbf{\gamma}}_n(s) \mathrm{d}s \,,
%> @f]
%> where the mapping @f$\hat{\mathbf{\gamma}}_n@f$ is given in 
%> <code>gammaMap()</code>.
%>
%> @param  N    The local number of degrees of freedom
%> @param  basesOnGamma  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi1D.
%> @retval ret  The computed array @f$[N\times N\times 3]@f$
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
function ret = integrateRefEdgePhiPhiHybrid(N, basesOnGamma)
validateattributes(basesOnGamma, {'struct'}, {}, mfilename, 'basesOnGamma')
p = N-1;  qOrd = 2*p+1;  [~, W] = quadRule1D(qOrd);
ret = zeros(N, N, 3); % [N x N x 3]
  for i = 1 : N
    for j = 1 : N
      ret(i, j) = sum( W' .* basesOnGamma.phi1D{qOrd}(:,i) .* basesOnGamma.phi1D{qOrd}(:,j) );
    end % for
  end % for
end % function
