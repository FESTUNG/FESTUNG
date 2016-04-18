% Evaluates all permutations of two basis functions
% in all quadrature points of all edges of the reference triangle and 
% multiplies with the according quadrature weight.
%
%===============================================================================
%> @file integrateRefEdgePhiIntPhiIntPerQuad.m
%>
%> @brief Evaluates all permutations of two basis functions
%>        in all quadrature points of all edges of the reference triangle
%>        and multiplies with the according quadrature weight.
%===============================================================================
%>
%> @brief Evaluates all permutations of two basis functions
%>        in all quadrature points of all edges of the reference triangle
%>        and multiplies with the according quadrature weight.
%>
%> It computes a multidimensional array @f$\hat{\mathsf{R}}^\mathrm{diag} 
%>    \in \mathbb{R}^\mathrm{diag}^{N\times N\times 3\times R}@f$, which is 
%> defined by
%> @f[
%> \left[\hat{\mathsf{R}}^\mathrm{diag}\right]_{i,j,n,r} \;:=\;
%> \hat{\varphi}_i\circ\hat{\mathbf{\gamma}}_n(q_r)\,
%> \hat{\varphi}_j\circ\hat{\mathbf{\gamma}}_n(q_r)\,
%> w_r
%> @f]
%> with the mapping @f$\hat{\mathbf{\gamma}}_n@f$ defined in 
%> <code>gammaMap()</code> and the quadrature points @f$q_r@f$ given by
%> <code>quadRule1D()</code>
%>
%> @param  N    The local number of degrees of freedom
%> @retval ret  The computed array @f$[N\times N\times 3\times R]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Florian Frank, Balthasar Reuter, Vadym Aizinger
%> Modified by Hennes Hajduk, 2016-04-18
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
function ret = integrateRefEdgePhiIntPhiIntPerQuad(N)
global gPhi1D
p = (sqrt(8*N+1)-3)/2;  qOrd = max(2*p+1,1);  [~, W] = quadRule1D(qOrd);
ret = zeros(N, N, 3, length(W)); % [N x N x 3 x R]
for n = 1 : 3 % 3 edges
  for i = 1 : N
    for j = 1 : N
      ret(i,j,n,:) = gPhi1D{qOrd}(:,i,n).* gPhi1D{qOrd}(:,j,n) .* W.';
    end % for
  end % for
end % for
end % function
