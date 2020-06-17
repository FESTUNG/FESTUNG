% Evaluates all permutations of two basis functions from different elements
% in all quadrature points of all edges of the reference triangle.

%===============================================================================
%> @file
%>
%> @brief Evaluates all permutations of two basis functions from different elements
%>        in all quadrature points of all edges of the reference triangle.
%===============================================================================
%>
%> @brief Evaluates all permutations of two basis functions from different elements
%>        in all quadrature points of all edges of the reference triangle.
%>
%> It computes a multidimensional array @f$\hat{\mathsf{R}}^\mathrm{offdiag} 
%>    \in \mathbb{R}^{N\times N\times 3\times 3\times R}@f$, which is defined by
%> @f[
%> \left[\hat{\mathsf{R}}^\mathrm{offdiag}\right]_{i,j,n^-,n^+,r} \;:=\;
%> \hat{\varphi}_i\circ\hat{\mathbf{\gamma}}_{n^-}(q_r)\,
%> \hat{\varphi}_j\circ\vartheta_{n^-n^+}\circ\hat{\mathbf{\gamma}}_{n^-}(q_r)
%> @f]
%> with the mapping @f$\hat{\mathbf{\gamma}}_n@f$ defined in 
%> <code>gammaMap()</code>, the mapping 
%> @f$\hat{\mathbf{\vartheta}}_{n^-n^+}@f$ as described in
%> <code>theta()</code>, and the quadrature points @f$q_r@f$ given by
%> <code>quadRule1D()</code>
%>
%> @param  N    The local number of degrees of freedom
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi1D and thetaPhi1D.
%> @param  qOrd (optional) The order of the quadrature rule to be used.
%> @retval ret  The computed array @f$[N\times N\times 3\times 3\times R]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2019 Florian Frank, Balthasar Reuter, Vadym Aizinger
%> @author Balthasar Reuter, 2019
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
function ret = integrateRefEdgeQuadriPhiIntPhiExtPerQuad(N, basesOnQuad, qOrd)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
[~, W] = quadRule1D(qOrd);
mapE0E = [2 1 4 3];
ret = zeros(N,N,4,length(W)); % [N x N x N x 4]
for nn = 1 : 4 % 4 edges
  np = mapE0E(nn);
  for i = 1 : N
    for j = 1 : N
      ret(i, j, nn, :) = basesOnQuad.phi1D{qOrd}(:,i,nn) .* basesOnQuad.phi1D{qOrd}(:,j,np) .* W.';
    end % for
  end % for
end % for
end
