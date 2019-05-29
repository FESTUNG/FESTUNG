% Compute integrals on the reference triangle, whose integrands consist of all
% permutations of a basis function with one linear Lagrange basis function.
%
%===============================================================================
%> @file integrateRefElemPhiLagrPhi.m
%>
%> @brief Compute integrals on the reference triangle, whose integrands 
%>        consist of all permutations of a basis function with one linear 
%>        Lagrange basis function.
%===============================================================================
%>
%> @brief Compute integrals on the reference triangle, whose integrands 
%>        consist of all permutations of a basis function with one linear 
%>        Lagrange basis function.
%>
%> It computes a multidimensional array
%> @f$\hat{\mathsf{M}_L} \in \mathbb{R}^{3 \times N}@f$
%> defined by
%> @f[
%>  [\hat{\mathsf{H}_L}]_{i,j} = \int_{\hat{T}} 
%>    \hat{\varphi}_i^L \hat{\varphi}_j \mathrm{d} \hat{\mathbf{x}}\,.
%> @f]
%> \hat{\varphi}_i^L is the piecewise linear function whose value in the 
%> i-th local vertex is one and zero in all others and \hat{\varphi}_j is
%> the j-th local DG basis function.
%>
%> @param  N    The local number of degrees of freedom
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. 
%> @retval ret  The computed array @f$[3\times N]@f$
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
function ret = integrateRefElemPhiLagrPhi(N, basesOnQuad)
ret = zeros(3, N);
p = (sqrt(8*N+1)-3)/2;  qOrd = max(2*p, 1);  [Q1, Q2, W] = quadRule2D(qOrd);
for i = 1 : 3
  switch i
    case 1
      phiLagr = 1 - Q1 - Q2;
    case 2
      phiLagr = Q1;
    case 3
      phiLagr = Q2;
    otherwise
      error('Invalid local index for Lagrange basis function.');
  end % switch
	for j = 1 : N
  	ret(i, j) = sum( phiLagr' .* W' .* basesOnQuad.phi2D{qOrd}(:,j) );
	end % for
end % for
end % function
