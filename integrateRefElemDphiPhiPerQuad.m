% Compute integral contributions on the reference triangle, whose integrands consist of all
% permutations of a basis function with one of the (spatial) derivatives of a
% basis function for every integration point.
%
%===============================================================================
%> @file integrateRefElemDphiPhiPerQuad.m
%>
%> @brief NEW Compute integral contributions on the reference triangle, whose integrands consist 
%>        of all permutations of a basis function with one of the (spatial)
%>        derivatives of a basis function  for every integration point.
%===============================================================================
%>
%> @brief Compute integral contributions on the reference triangle @f$\hat{T}@f$,
%>        whose integrands consist of all permutations of a basis function with
%>        one of the (spatial) derivatives of a basis function for every integration point @f$ {\boldsymbol{\hat{q}}}_r @f$.
%>
%> It computes a multidimensional array
%> @f$\hat{\mathsf{G}} \in \mathbb{R}^{2 \times N \times N \times R}@f$
%> defined by
%> @f[
%> 	[\mathsf{\hat{G}}]_{m,i,j,r} = \omega_{r} \, \partial_{\hat{x}^{m}}{} \, \hat{\varphi}_{ki}({\boldsymbol{\hat{q}}}_r) \, \hat{\varphi}_{kj}({\boldsymbol{\hat{q}}}_r)
%> @f]
%> with @f$\omega_{r}@f$ denoting the integration weight. For an efficient access of the array entris it is actually stored as an @f$2 \times 1 \text{ cell}@f$ with each cell storing a @f$N \times N \times R@f$ array.
%> 
%> 
%> @param  N    The local number of degrees of freedom @f$[\text{scalar}]@f$
%> 
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi2D and gradPhi2D. @f$[\text{scalar}]@f$
%> 
%> @param  qord 	(optional) The order of the quadrature rule to be used.. @f$[\text{scalar}]@f$
%> 
%> @retval ret  The computed array @f$2 \times 1 \text{ cell}@f$ 
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
%> 
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
%%TODO 
%
% I precompute matrices G_bar that allows for an 'easy' evaluation of 
% u_{m} phi_{kj} \partial_{x_m} phi_{ki}
%
function ret = integrateRefElemDphiPhiPerQuad(N, basesOnQuad, qOrd)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
if nargin < 3, p = (sqrt(8*N+1)-3)/2;  qOrd = max(2*p, 1);  end
[~, ~, W] = quadRule2D(qOrd); R = length(W);
ret = { zeros(N, N, R); zeros(N, N, R) };
for i = 1 : N
  for j = 1 : N
    for m = 1 : 2
      ret{m}(i, j, :) =  W(:) .* basesOnQuad.phi2D{qOrd}(:, j) .* basesOnQuad.gradPhi2D{qOrd}(:, i, m);
    end % for
  end % for
end % for
end % function
