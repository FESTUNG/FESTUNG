% Compute integrals over the edges of the reference triangle, whose integrands 
% consist of all permutations of two basis functions from different elements.
%
%===============================================================================
%> @file integrateRefEdgePhiIntPhiExt.m
%>
%> @brief Compute integrals over the edges of the reference triangle, whose 
%>        integrands consist of all permutations of two basis functions from 
%>        different elements.
%===============================================================================
%>
%> @brief Compute integrals over the edges of the reference triangle 
%>        @f$\hat{T}@f$, whose integrands consist of all permutations of two
%>        basis functions, of which one belongs to a neighboring element that
%>        is transformed using @f$\hat{\mathbf{\vartheta}}@f$.
%>
%> It computes a multidimensional array @f$\hat{\mathsf{{S}}}^\mathrm{offdiag} 
%>    \in \mathbb{R}^{N\times N\times 3\times 3}@f$, which is defined by
%> @f[
%> [\hat{\mathsf{{S}}}^\mathrm{offdiag}]_{i,j,n^-,n^+} =
%>   \int_0^1 \hat{\varphi}_i \circ \hat{\mathbf{\gamma}}_{n^-}(s) 
%>   \hat{\varphi}_j\circ \hat{\mathbf{\vartheta}}_{n^-n^+} \circ
%>   \hat{\mathbf{\gamma}}_{n^-}(s) \mathrm{d}s \,,
%> @f]
%> with the mapping @f$\hat{\mathbf{\gamma}}_n@f$ defined in 
%> <code>gammaMap()</code> and the mapping 
%> @f$\hat{\mathbf{\vartheta}}_{n^-n^+}@f$ as described in <code>theta()</code>.
%>
%> @param  N    The local number of degrees of freedom
%> @retval ret  The computed array @f$[N\times N\times 3\times 3]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2015 Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function ret = integrateRefEdgePhiIntPhiExt(N)
global gPhi1D gThetaPhi1D
p = (sqrt(8*N+1)-3)/2;  qOrd = 2*p+1;  [~, W] = quadRule1D(qOrd); 
ret = zeros(N, N, 3, 3); % [N x N x 3 x 3]
for nn = 1 : 3
  for np = 1 : 3
    for i = 1 : N
      for j = 1 : N
        ret(i, j, nn, np) = sum( W' .* gPhi1D{qOrd}(:,i,nn) .* gThetaPhi1D{qOrd}(:,j,nn,np) );
      end % for
    end % for
  end % for
end % for
end % function
