% Evaluates the linear Lagrange basis functions in all quadrature points of
% all edges of the reference triangle.
%===============================================================================
%> @file integrateRefEdgePhiLagrPerQuad.m
%>
%> @brief Evaluates the linear Lagrange basis functions in all quadrature 
%>        points of all edges of the reference triangle.
%===============================================================================
%>
%> @brief Evaluates the linear Lagrange basis functions in all quadrature 
%>        points of all edges of the reference triangle.
%>
%> It computes a multidimensional array
%> @f$\hat{\mathsf{{S}_L}}}\in\mathbb{R}^{3\times 3\times R}@f$
%> defined by
%> @f[
%> [\hat{\mathsf{{S}_L}}]_{i,n,r} =
%>   \hat{omega}_r \hat{\varphi}_i^L \circ \hat{\mathbf{\gamma}}_n(\hat{q}_r)\,,
%> @f]
%> where the mapping @f$\hat{\mathbf{\gamma}}_n@f$ is given in 
%> <code>gammaMap()</code> and \hat{\varphi}_i^L is the piecewise linear 
%> function whose value in the i-th local vertex is one and zero in all 
%> others.
%>
%> @param  N    The local number of degrees of freedom
%> @retval ret  The computed array @f$[3\times 3\times R]@f$
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
function ret = integrateRefEdgePhiLagrPerQuad(N)
p = (sqrt(8*N+1)-3)/2;  qOrd = 2*p+1;  [Q, W] = quadRule1D(qOrd);
ret = zeros(3, 3, length(W)); % [3 x 3 x R]
for n = 1 : 3 % 3 edges
	[X, Y] = gammaMap(n, Q);
  ret(:,n,:) = [phiLagr(1,X,Y); phiLagr(2,X,Y); phiLagr(3,X,Y)] .* repmat(W, 3, 1);
end % for
end % function
