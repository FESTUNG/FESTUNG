% Compute integrals over the reference edge, whose integrands 
% consist of all permutations of two edge basis functions.

%===============================================================================
%> @file
%>
%> @brief Compute integrals over the reference edge, whose integrands 
%> consist of all permutations of two edge basis functions.
%===============================================================================
%>
%> @brief Compute integrals over the edges of the reference edge 
%>        @f$\hat{E}@f$, whose integrands consist of all permutations of two
%>        edge basis functions.
%>
%> It computes a matrix
%> @f$\hat{\mathsf{M}}_{\mu}\in\mathbb{R}^{\bar{N}\times \bar{N}}@f$
%> defined by
%> @f[
%> [\hat{\mathsf{M}}_{\mu}]_{i,j} =
%>   \int_{\hat{E}} \hat{\mu}_i \hat{\mu}_j \mathrm{d}\hat{\mathbf{s}} 
%>   =
%>   \int_{0}^{1} \hat{\mu}_i \hat{\mu}_j \mathrm{d}{\mathbf{s}} \,.\
%> @f]
%>
%> @param  N    The local number of degrees of freedom
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least mu. @f$[\text{struct}]@f$
%> @param  qOrd 	(optional) The order of the quadrature rule to be used.. @f$[\text{scalar}]@f$
%> @retval ret  The computed matrix @f$[\bar{N} \times \bar{N}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
%> @author Alexander Jaust, 2017
%> @author Balthasar Reuter, 2017
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
function ret = integrateRefEdgeMuMu(N, basesOnQuad, qOrd)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnGamma')
if nargin < 3, p = N-1;  qOrd = 2*p+1;  end
[~, W] = quadRule1D(qOrd);
ret = zeros(N, N); % [N x N]
for i = 1 : N
  for j = 1 : N
    ret(i, j) = sum( W' .* basesOnQuad.mu{qOrd}(:,i) .* basesOnQuad.mu{qOrd}(:,j) );
  end % for
end % for
end % function
