% Compute integrals over the edges of the reference triangle, whose integrands 
% consist of all permutations of two linear Lagrange basis functions.

%===============================================================================
%> @file integrateRefElemPhiLagrPhiLagr.m
%>
%> @brief Compute integrals over the edges of the reference triangle, whose 
%>        integrands consist of all permutations of two linear Lagrange basis 
%>        functions.
%===============================================================================
%>
%> @brief Compute integrals over the edges of the reference triangle, whose 
%>        integrands consist of all permutations of two linear Lagrange basis 
%>        functions.
%>
%> It computes a matrix
%> @f$\hat{\mathsf{{M}_L}}\in\mathbb{R}^{3\times 3}@f$
%> defined by
%> @f[
%> [\hat{\mathsf{{M}_L}}_{i,j} =
%>   \int_{\hat{T}} \hat{\varphi}_i^L \hat{\varphi}_j^L \mathrm{d} \hat{\mathbf{x}}\,.
%> @f]
%> where \hat{\varphi}_i^L is the piecewise linear function whose value in the 
%> i-th local vertex is one and zero in all others and analogously for \hat{\varphi}_j.
%>
%> @retval ret  The computed array @f$[3\times 3]@f$
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
function ret = integrateRefElemPhiLagrPhiLagr()
[Q1, Q2, W] = quadRule2D(2);
ret = zeros(3,3); % 3 x 3
for i = 1 : 3
  for j = 1 : 3
    ret(i,j) = sum(phiLagr(i,Q1,Q2) .* phiLagr(j,Q1,Q2) .* W);
  end % for
end % for
end % functiom
