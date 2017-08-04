% Evaluates all permutations of two basis functions
% in all quadrature points of all edges of the reference triangle.

%===============================================================================
%> @file integrateRefEdgePhiIntPhiIntPerQuad.m
%>
%> @brief Evaluates all permutations of two basis functions
%>        in all quadrature points of all edges of the reference triangle.
%===============================================================================
%>
%> @brief Evaluates all permutations of two basis functions
%>        in all quadrature points of all edges of the reference triangle.
%>
%> It computes a multidimensional array @f$\hat{\mathsf{R}}^\mathrm{diag} 
%>    \in \mathbb{R}^{N\times N\times 3\times R}@f$, which is defined by
%> @f[
%> \left[\hat{\mathsf{R}}^\mathrm{diag}\right]_{i,j,n,r} \;:=\;
%> \hat{\varphi}_i\circ\hat{\mathbf{\gamma}}_n(q_r)\,
%> \hat{\varphi}_j\circ\hat{\mathbf{\gamma}}_n(q_r)\,
%> @f]
%> with the mapping @f$\hat{\mathbf{\gamma}}_n@f$ defined in 
%> <code>gammaMap()</code> and the quadrature points @f$q_r@f$ given by
%> <code>quadRule1D()</code>
%>
%> @param  N    The local number of degrees of freedom
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi1D.
%> @param  qOrd (optional) The order of the quadrature rule to be used.
%> @retval ret  The computed array @f$[N\times N\times 3\times R]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function ret = integrateRefEdgePhiIntPhiIntPerQuad(N, basesOnQuad, qOrd)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
if nargin < 3
  p = (sqrt(8*N+1)-3)/2;  qOrd = max(2*p+1,1);
end % if
[~, W] = quadRule1D(qOrd);
ret = zeros(N, N, 3, length(W)); % [N x N x 3 x R]
for n = 1 : 3 % 3 edges
  for i = 1 : N
    for j = 1 : N
      ret(i,j,n,:) = basesOnQuad.phi1D{qOrd}(:,i,n) .* basesOnQuad.phi1D{qOrd}(:,j,n) .* W.';
    end % for
  end % for
end % for
end % function
