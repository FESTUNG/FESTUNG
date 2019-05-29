% Compute integrals on the reference triangle, whose integrands consist of all
% permutations of two (spatial) derivatives of linear Lagrange basis functions.
%
%===============================================================================
%> @file integrateRefElemDphiLagrDphiLagr.m
%>
%> @brief Compute integrals on the reference triangle, whose integrands consist
%>        of all permutations of two (spatial) derivatives of linear Lagrange 
%>        basis functions.
%===============================================================================
%>
%> @brief Compute integrals on the reference triangle, whose integrands consist
%>        of all permutations of two (spatial) derivatives of linear Lagrange 
%>        basis functions.
%>
%> It computes a multidimensional array
%> @f$\hat{\mathsf{A}_L} \in \mathbb{R}^{3 \times 3 \times 2 \times 2}@f$
%> defined by
%> @f[
%>  [\hat{\mathsf{H}_L}]_{i,j,m_1,m_2} = \int_{\hat{T}} \partial_{\hat{x}^{m_1}} 
%>    \hat{\varphi}_i^L \partial_{\hat{x}^{m_2}} \hat{\varphi}_j^L 
%>    \mathrm{d} \hat{\mathbf{x}}\,.
%> @f]
%> \hat{\varphi}_{i,j}^L is the piecewise linear function whose value in the 
%> i/j-th local vertex is one and zero in all others.
%>
%> @param  N    The local number of degrees of freedom
%> @retval ret  The computed array @f$[3\times 3\times 2\times 2]@f$
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
function ret = integrateRefElemDphiLagrDphiLagr(N)
ret = zeros(3, 3, 2, 2);
p = (sqrt(8*N+1)-3)/2;  qOrd = max(2*p,1);  [Q1, Q2, W] = quadRule2D(qOrd);
for i = 1 : 3
	for j = 1 : 3
		for m1 = 1 : 2
      for m2 = 1 : 2
  			ret(i, j, m1, m2) = sum( gradPhiLagr(i,m1,Q1,Q2) .* gradPhiLagr(j,m2,Q1,Q2) .* W );
      end % for
		end % for
	end % for
end % for
end % function
