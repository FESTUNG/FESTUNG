% Compute integrals over the edges of the reference triangle, whose integrands 
% consist of one basis function.

%===============================================================================
%> @file
%>
%> @brief Compute integrals over the edges of the reference triangle, whose 
%>        integrands consist of one basis function.
%===============================================================================
%>
%> @brief Compute integrals over the edges of the reference triangle 
%>        @f$\hat{T}@f$, whose integrands consist of one basis function.
%>
%> It computes a matrix @f$\hat{\mathsf{{R}}}
%> \in \mathbb{R}^{N\times 3}@f$, which is defined by
%> @f[
%> [\hat{\mathsf{{R}}}]_{j,n} =
%>   \int_0^1 \hat{\varphi}_j\circ \hat{\mathbf{\gamma}}_{n^-}(s)\mathrm{d}s \,,
%> @f]
%> with the mapping @f$\hat{\mathbf{\gamma}}_n@f$ defined in 
%> <code>gammaMap()</code>.
%>
%> @param  N            The local number of degrees of freedom.
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi1D.
%> @retval ret  The computed array @f$[N\times 3]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2018 Florian Frank, Balthasar Reuter, Vadym Aizinger
%> author Hennes Hajduk, 2018
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
function ret = integrateRefEdgePhiInt(N, basesOnQuad)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')

ret = zeros(N, 3);
p = (sqrt(8*N+1)-3)/2;  qOrd = 2*p+1;  [~, W] = quadRule1D(qOrd);
for nn = 1 : 3 % 3 edges
  for j = 1 : N
    ret(j, nn) = sum( W'.* basesOnQuad.phi1D{qOrd}(:,j,nn) );
  end % for
end % for
end % function
