% Evaluate all permutations of three basis functions in the vertices of the 
% unit interval.

%===============================================================================
%> @file computePhiIntPhiIntPhiIntV0T1D.m
%>
%> @brief Evaluate all permutations of three basis functions in the vertices 
%>        of the unit interval.
%===============================================================================
%>
%> @brief Evaluate all permutations of three basis functions in the vertices 
%>        of the unit interval.
%>
%> It computes two multidimensional array @f$\hat{\mathsf{{R}}}^s
%> \in \mathbb{R}^{N\times N\times N\times 2}, s\in\{1,2\}@f$, which is defined by
%> @f[
%> [\hat{\mathsf{{R}}}^1]_{i,j,l,n} =
%>   \hat{\phi}_i (n-1) \hat{\phi}_l (n-1) \hat{\phi}_j (n-1) \,,
%> [\hat{\mathsf{{R}}}^2]_{i,j,l,n} =
%>   \hat{\phi}_i (n-1) \hat{\phi}_l (n-1) \hat{\phi}_j (n-1) (n-1) \,.
%> @f]
%>
%> @param  N            The local number of degrees of freedom
%> @param  qOrd         The order of the quadrature rule to be used.
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi0D.
%> @retval ret  The computed arrays @f$[2\times1 \text{cell}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
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
function ret = computePhiIntPhiIntPhiIntV0T1D(N, qOrd, basesOnQuad)
ret = { zeros(N, N, N, 2), zeros(N, N, N, 2) };
for n = 1 : 2
  for i = 1 : N
    for j = 1 : N
      for l = 1 : N
        ret{1}(i,j,l,n) = basesOnQuad.phi0D{qOrd}(i,n) * basesOnQuad.phi0D{qOrd}(j,n) * basesOnQuad.phi0D{qOrd}(l,n);
        ret{2}(i,j,l,n) = (n-1) * basesOnQuad.phi0D{qOrd}(i,n) * basesOnQuad.phi0D{qOrd}(j,n) * basesOnQuad.phi0D{qOrd}(l,n);
      end % for l
    end % for j
  end  % for i
end  % for n
end  % function