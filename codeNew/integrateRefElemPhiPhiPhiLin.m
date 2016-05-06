% Compute integrals on the reference triangle, whose integrands 
% consist of all permutations of two basis functions and a linear basis function.
%
%===============================================================================
%> @file integrateRefElemPhiPhiPhiLin.m
%>
%> @brief Compute integrals on the reference triangle, whose 
%>        integrands consist of all permutations of two basis functions and a 
%>        linear basis function.
%===============================================================================
%>
%> @brief Compute integrals on the reference triangle @f$\hat{T}@f$,
%>        whose integrands consist of all permutations of two basis functions 
%>        and a linear basis function.
%>
%> It computes a multidimensional array
%> @f$\hat{\mathsf{D}}\in\mathbb{R}^{N\times N\times 3}@f$
%> defined by
%> @f[
%> [\hat{\mathsf{M}}]_{i,j,l} =
%>   \int_{\hat{T}} \hat{\varphi}_i \hat{\varphi}_j hat{\varphi}_l \mathrm{d}\hat{\mathbf{x}} \,.
%> @f]
%>
%> @param  N    The local number of degrees of freedom
%> @retval ret  The computed matrix @f$[N\times N \times 3]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Hennes Hajduk, Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function ret = integrateRefElemPhiPhiPhiLin(N)
global gPhi2D
p = (sqrt(8*N+1)-3)/2; qOrd = max(2*p, 1); [Q1, Q2, W] = quadRule2D(qOrd);
ret = zeros(N,N,3);
for i = 1 : N
  for j = 1 : N
    for l = 1 : 3
      ret(i,j,l) = sum( W.' .* gPhi2D{qOrd}(:,i) .* gPhi2D{qOrd}(:,j) .* phi(l, Q1, Q2).' );
    end % for
  end % for
end % for
end % function