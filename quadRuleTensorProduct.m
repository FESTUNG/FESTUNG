% Provides quadrature points and associated weights built from tensor
% product of one-dimensional quadrature rules.

%===============================================================================
%> @file quadRuleTensorProduct.m
%>
%> @brief Provides quadrature points and associated weights built from tensor
%>        product of one-dimensional quadrature rules.
%===============================================================================
%>
%> @brief Provides quadrature points and associated weights built from tensor
%>        product of one-dimensional quadrature rules.
%>
%> <code>[Q1, Q2, W] = quadRuleTensorProduct(qOrd)</code> returns quadrature points 
%> @f$\hat{\vec{q}}_r = [\hat{q}_r^1, \hat{q}_r^2]^\mathrm{T}@f$ within the reference
%> square @f$\hat{T}@f$ in lists of @f$\hat{x}^1@f$ and @f$\hat{x}^2@f$ coordinates
%> <code>Q1</code> and <code>Q2</code>, respectively, and the associated 
%> weights <code>W</code>.  The quadrature rule is exact for polynomials of order 
%> <code>qOrd</code> (see @ref FRAK2015 for details).
%>
%> If no one-dimensional quadrature-rules are specified, they default to 
%> <code>quadRule1D()</code>.
%> 
%> @param  qOrd The order of the quadrature rule.
%> @param  quadRule1  (optional) Function handle for one-dimensional 
%>                    quadrature-rule in @f$\hat{x}^1@f$-direction.
%> @param  quadRule2  (optional) Function handle for one-dimensional 
%>                    quadrature-rule in @f$\hat{x}^2@f$-direction.
%> @retval Q1   The @f$\hat{x}^1@f$ coordinates of the quadrature points.
%> @retval Q2   The @f$\hat{x}^2@f$ coordinates of the quadrature points.
%> @retval W    The associated weights.
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
function [Q1, Q2, W] = quadRuleTensorProduct(qOrd, quadRule1, quadRule2)
if nargin < 2
  quadRule1 = @quadRule1D;
  quadRule2 = @quadRule1D;
end % if
[Q1, W1] = quadRule1(qOrd);
[Q2, W2] = quadRule2(qOrd);
[Q1, Q2] = meshgrid(Q1, Q2);
[W1, W2] = meshgrid(W1, W2);
Q1 = Q1(:)'; Q2 = Q2(:)'; W = reshape(W1 .* W2, [], 1)';
end % function