% Compute integrals on the reference triangle, whose integrands consist of all
% permutations of a basis function with one of the (spatial) derivatives of a
% linear Lagrange basis function.
%
%===============================================================================
%> @file integrateRefElemDphiLagrPhi.m
%>
%> @brief Compute integrals on the reference triangle, whose integrands consist 
%>				of all permutations of a basis function with one of the (spatial) 
%>				derivatives of a linear Lagrange basis function.
%===============================================================================
%>
%> @brief Compute integrals on the reference triangle, whose integrands consist 
%>				of all permutations of a basis function with one of the (spatial) 
%>				derivatives of a linear Lagrange basis function.
%>
%> It computes a multidimensional array
%> @f$\hat{\mathsf{H}_L} \in \mathbb{R}^{3 \times N \times 2}@f$
%> defined by
%> @f[
%>  [\hat{\mathsf{H}_L}]_{i,j,m} = \int_{\hat{T}} \partial_{\hat{x}^m} 
%>    \hat{\varphi}_i^L \hat{\varphi}_j \mathrm{d} \hat{\mathbf{x}}\,.
%> @f]
%> \hat{\varphi}_i^L is the piecewise linear function whose value in the 
%> i-th local vertex is one and zero in all others and \hat{\varphi}_j is
%> the j-th local DG basis function.
%>
%> @param  N    The local number of degrees of freedom
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points.
%> @retval ret  The computed array @f$[3\times N\times 2]@f$
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
function ret = integrateRefElemDphiLagrPhi(N, basesOnQuad)
ret = zeros(3, N, 2);
p = (sqrt(8*N+1)-3)/2;  qOrd = max(2*p, 1);  [Q1,Q2,W] = quadRule2D(qOrd);
for i = 1 : 3
	for j = 1 : N
		for m = 1 : 2
			ret(i, j, m) = sum( gradPhiLagr(i,m,Q1,Q2)' .* W' .* basesOnQuad.phi2D{qOrd}(:,j) );
		end % for
	end % for
end % for
end % function
