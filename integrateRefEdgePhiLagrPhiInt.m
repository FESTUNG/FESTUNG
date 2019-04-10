% Compute integrals over the edges of the reference triangle, whose integrands 
% consist of all permutations of a linear Lagrange basis function and a DG 
% basis function.

%===============================================================================
%> @file integrateRefEdgePhiLagrPhiInt.m
%>
%> @brief Compute integrals over the edges of the reference triangle, whose
%>        integrands consist of all permutations of a linear Lagrange basis
%>        function and a DG basis function.
%===============================================================================
%>
%> @brief Compute integrals over the edges of the reference triangle, whose
%>        integrands consist of all permutations of a linear Lagrange basis
%>        function and a DG basis function.
%>
%> It computes a multidimensional array
%> @f$\hat{\mathsf{{S}_L}}^\mathrm{diag}\in\mathbb{R}^{3\times N\times3}@f$
%> defined by
%> @f[
%> [\hat{\mathsf{{S}_L}}^\mathrm{diag}]_{i,j,n} =
%>   \int_0^1 \hat{\varphi}_i^L \circ \hat{\mathbf{\gamma}}_n(s) 
%>   \hat{\varphi}_j\circ \hat{\mathbf{\gamma}}_n(s) \mathrm{d}s \,,
%> @f]
%> where the mapping @f$\hat{\mathbf{\gamma}}_n@f$ is given in 
%> <code>gammaMap()</code> and \hat{\varphi}_i^L is the piecewise linear 
%> function whose value in the i-th local vertex is one and zero in all 
%> others and  \hat{\varphi}_j is the j-th local DG basis function.
%>
%> @param  N    The local number of degrees of freedom
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi1D.
%> @retval ret  The computed array @f$[N\times N\times 3]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Hennes Hajduk, Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function ret = integrateRefEdgePhiLagrPhiInt(N, basesOnQuad)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
p = (sqrt(8*N+1)-3)/2;  qOrd = 2*p+1;  [Q, W] = quadRule1D(qOrd);
ret = zeros(3, N, 3); % [3 x N x 3]
for n = 1 : 3 % 3 edges
	[X, Y] = gammaMap(n, Q);
  for i = 1 : 3
    for j = 1 : N
      ret(i, j, n) = sum( W' .* phiLagr(i, X, Y)' .* basesOnQuad.phi1D{qOrd}(:,j,n) );
    end % for
  end % for
end % for
end % function
