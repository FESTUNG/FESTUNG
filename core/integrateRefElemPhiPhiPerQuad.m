% Compute integral contributions on the reference triangle, whose integrands consist of all
% permutations of two basis functions for every integration point.

%===============================================================================
%> @file
%>
%> @brief Compute integral contributions on the reference triangle, whose integrands consist 
%>        of all permutations of two basis functions for every integration point.
%===============================================================================
%>
%> @brief Compute integral contributions on the reference triangle @f$\hat{T}@f$,
%>        whose integrands consist of all permutations of two basis functions
%>        for every integration point @f$ {\boldsymbol{\hat{q}}}_r @f$.
%>
%> It computes a multidimensional array
%> @f$\hat{\mathsf{G}} \in \mathbb{R}^{N \times N \times R}@f$
%> defined by
%> @f[
%> 	[\mathsf{\hat{G}}]_{i,j,r} = \omega_{r} \, \hat{\varphi}_{ki}({\boldsymbol{\hat{q}}}_r) \, \hat{\varphi}_{kj}({\boldsymbol{\hat{q}}}_r)
%> @f]
%> with @f$\omega_{r}@f$ denoting the integration weight.
%> 
%>
%> @param  N    The local number of degrees of freedom @f$[\text{scalar}]@f$
%> 
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi2D. @f$[\text{struct}]@f$
%> 
%> @param  qOrd 	(optional) The order of the quadrature rule to be used. @f$[\text{scalar}]@f$
%> 
%> @retval ret  The computed array @f$[N \times N \times R]@f$ 
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2019 Balthasar Reuter, Florian Frank, Vadym Aizinger
%> 
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
function ret = integrateRefElemPhiPhiPerQuad(N, basesOnQuad, qOrd)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
if nargin < 3, p = (sqrt(8*N+1)-3)/2;  qOrd = max(2*p, 1);  end
[~, ~, W] = quadRule2D(qOrd); R = length(W);
ret = zeros(N, N, R);
for i = 1 : N
  for j = 1 : N
    ret(i, j, :) =  W(:) .* basesOnQuad.phi2D{qOrd}(:, i) .* basesOnQuad.phi2D{qOrd}(:, j);
  end % for
end % for
end % function
