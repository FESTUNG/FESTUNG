% Compute integrals over the edges of the reference square, whose integrands 
% consist of all permutations of two basis functions.

%===============================================================================
%> @file integrateRefEdgeTetraPhiIntPhiInt.m
%>
%> @brief Compute integrals over the edges of the reference square, whose 
%>        integrands consist of all permutations of two basis functions.
%===============================================================================
%>
%> @brief Compute integrals over the edges of the reference square 
%>        @f$\hat{T}@f$, whose integrands consist of all permutations of two
%>        basis functions.
%>
%> It computes a multidimensional array @f$\hat{\mathsf{{S}}}^\mathrm{offdiag} 
%>    \in \mathbb{R}^{N\times N\times 4}@f$, which is defined by
%> @f[
%> [\hat{\mathsf{{S}}}^\mathrm{offdiag}]_{i,j,n^-} =
%>   \int_0^1 \hat{\varphi}_i \circ \hat{\mathbf{\gamma}}_{n^-}(s) 
%>   \hat{\varphi}_j\circ \hat{\mathbf{\gamma}}_{n^-}(s) \mathrm{d}s \,,
%> @f]
%> with the mapping @f$\hat{\mathbf{\gamma}}_n@f$ defined in 
%> <code>gammaMapTetra()</code>.
%>
%> @param  N            The local number of degrees of freedom.
%> @param  qOrd         The order of the quadrature rule to be used.
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi1D.
%> @retval ret  The computed array @f$[N\times N\times 4]@f$
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
function ret = integrateRefEdgeTetraPhiIntPhiInt(N, qOrd, basesOnQuad)
[~, W] = quadRule1D(qOrd);
ret = zeros(N, N, 4);
for n = 1 : 4
  for i = 1 : N
    for j = 1 : i
      ind = sub2ind([N N 4], [i j], [j i], [n n]);
      ret(ind) = W * ( basesOnQuad.phi1D{qOrd}(:,i,n) .* basesOnQuad.phi1D{qOrd}(:,j,n) );
    end  % for j
  end  % for i
end  % for n
end  % function