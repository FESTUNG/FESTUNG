% Provides quadrature points and associated weights within the reference triangle.
%
%===============================================================================
%> @file quadRule2D.m
%>
%> @brief Provides quadrature points and associated weights within the reference 
%>        triangle.
%===============================================================================
%>
%> @brief Provides quadrature points and associated weights within the reference 
%>        triangle.
%>
%> Quadrature rules up to order 6 are implemented directly. For all higher orders
%> we call <code>triquad()</code> that uses Gaussian quadrature points on a
%> square which is collapsed to a triangle.  The area @f$1/2@f$ of the reference
%> triangle @f$\hat{T}@f$ is incorporated in the weights.
%>
%> <code>[Q1, Q2, W] = quadRule2D(qOrd)</code> returns quadrature points 
%> @f$\hat{\vec{q}}_r = [\hat{q}_r^1, \hat{q}_r^2]^\mathrm{T}@f$ within the reference triangle @f$\hat{T}@f$ in lists of
%> @f$\hat{x}^1@f$ and @f$\hat{x}^2@f$ coordinates <code>Q1</code> and <code>Q2</code>, respectively, and the associated 
%> weights <code>Q</code>.  The quadrature rule is exact for polynomials of order <code>qOrd</code> (see @ref FRAK2015 for details).
%> The area of @f$1/2@f$ of the reference triangle @f$\hat{T}@f$ is incorporated in the weights.
%>
%> An overview of quadrature rules on triangles is found in the &ldquo;Encyclopaedia of Cubature Formulas&rdquo; @ref Cools2003.
%>
%> @par Example
%> @code
%> [Q1, Q2, W] = quadRule2D(2);
%> N = 3; M = eye(N);
%> for i = 1:N
%>   for j = 1:N
%>     M(i,j) = sum(W .* phi(i,Q1,Q2) .* phi(j,Q1,Q2));
%>   end
%> end
%> @endcode
%>
%> 
%> 
%> @param  qOrd The order of the quadrature rule.
%> @retval Q1   The @f$\hat{x}^1@f$ coordinates of the quadrature points.
%> @retval Q2   The @f$\hat{x}^2@f$ coordinates of the quadrature points.
%> @retval W    The associated weights.
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
function [Q1, Q2, W] = quadRule2D(qOrd)
switch qOrd
  case {0, 1} % R = 1
    Q1 = 1/3;              Q2 = 1/3;              W  = 1/2;
  case 2 % R = 3
    Q1 = [1/6, 2/3, 1/6];  Q2 = [1/6, 1/6, 2/3];  W  = [1/6, 1/6, 1/6];
  case 3 % R = 4
    Q1 = [0.666390246, 0.178558728, 0.280019915, 0.075031109];
    Q2 = [0.178558728, 0.666390246, 0.075031109, 0.280019915];
    W  = [0.159020691, 0.159020691, 0.090979309, 0.090979309];
  case 4 % R = 6
    Q1 = [0.445948490915965, 0.108103018168070, 0.445948490915965, ...
          0.091576213509771, 0.816847572980458, 0.091576213509771];
    Q2 = [0.108103018168070, 0.445948490915965, 0.445948490915965, ...
          0.816847572980458, 0.091576213509771, 0.091576213509771];
    W  = [0.111690794839005, 0.111690794839005, 0.111690794839005, ...
          0.054975871827661, 0.054975871827661, 0.054975871827661];
  case 5 % R = 7
    a1 = (6-sqrt(15))/21;     a2 = (6+sqrt(15))/21;
    w1 = (155-sqrt(15))/2400; w2 = (155+sqrt(15))/2400;
    Q1 = [1/3,     a1, 1-2*a1,     a1,     a2, 1-2*a2,     a2];
    Q2 = [1/3, 1-2*a1,     a1,     a1, 1-2*a2,     a2,     a2];
    W  = [9/80,    w1,     w1,     w1,     w2,     w2,     w2];
  case 6 % R = 12
    Q1 = [0.063089014491502, 0.873821971016996, 0.063089014491502, ...
          0.249286745170910, 0.501426509658179, 0.249286745170910, ...
          0.310352451033785, 0.053145049844816, 0.636502499121399, ...
          0.053145049844816, 0.636502499121399, 0.310352451033785];
    Q2 = [0.063089014491502, 0.063089014491502, 0.873821971016996, ...
          0.249286745170910, 0.249286745170910, 0.501426509658179, ...
          0.053145049844816, 0.310352451033785, 0.053145049844816, ...
          0.636502499121399, 0.310352451033785, 0.636502499121399];
    W  = [0.025422453185103, 0.025422453185103, 0.025422453185103, ...
          0.058393137863189, 0.058393137863189, 0.058393137863189, ...
          0.041425537809187, 0.041425537809187, 0.041425537809187, ...
          0.041425537809187, 0.041425537809187, 0.041425537809187];
  otherwise % use Gauss quadrature points on square
    [X,Y,Wx,Wy] = triquad(qOrd, [0 0; 0 1; 1 0]); % third party function, see references
    Q1 = X(:).';
    Q2 = Y(:).';
    W = Wx * Wy.';  W = W(:).';
end % switch
end % function
