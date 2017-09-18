% Provides Gauss-Legendre quadrature points and associated weights on (0,1).

%===============================================================================
%> @file quadRule1D.m
%>
%> @brief Provides Gauss-Legendre quadrature points and associated weights on (0,1).
%===============================================================================
%>
%> @brief Provides a list of Gauss-Legendre quadrature points within the
%>        interval @f$[0,1]@f$ and associated positive weights.
%>
%> The quadrature rule is exact for polynomials up to the given order. The
%> length of the interval @f$(0,1)@f$ is incorporated in the weights. 
%>
%> @note A rule with @f$n@f$ points is exact for polynomials up to order @f$2n-1@f$.
%>
%> @par Example
%> @code
%> f = @(s) s.^2;
%> [Q, W] = quadRule1D(2);
%> F = dot(W, f(Q));
%> @endcode
%> 
%> @param  qOrd The order of the quadrature rule.
%> @retval Q    The list of Gauss-Legendre quadrature points within @f$(0,1)@f$.
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
function [Q, W] = quadRule1D(qOrd)
switch qOrd
  case {0, 1} % R = 1, number of quadrature points
    Q = 0;
    W = 2;
  case {2, 3} % R = 2
    Q = sqrt(1/3)*[-1, 1];
    W = [1, 1];
  case {4, 5} % R = 3
    Q = sqrt(3/5)*[-1, 0, 1];
    W = 1/9*[5, 8, 5];
  case {6, 7} % R = 4
    Q = [-1,-1,1,1].*sqrt(3/7+[1,-1,-1,1]*2/7*sqrt(6/5));
    W = 1/36*(18 + sqrt(30)*[-1,1,1,-1]);
  case {8, 9} % R = 5
    Q = [-1,-1,0,1,1].*sqrt(5+[2,-2,0,-2,2]*sqrt(10/7))/3;
    W = 1/900*(322+13*sqrt(70)*[-1,1,0,1,-1]+[0,0,190,0,0]);
  case {10, 11} % R = 6
    Q = [ 0.6612093864662645, -0.6612093864662645, -0.2386191860831969, ...
          0.2386191860831969, -0.9324695142031521,  0.9324695142031521];
    W = [ 0.3607615730481386,  0.3607615730481386,  0.4679139345726910, ...
          0.4679139345726910,  0.1713244923791704,  0.171324492379170];
  case {12, 13} % R = 7
    Q = [ 0.0000000000000000,  0.4058451513773972, -0.4058451513773972, ...
         -0.7415311855993945,  0.7415311855993945, -0.9491079123427585, ...
          0.9491079123427585];
    W = [ 0.4179591836734694,  0.3818300505051189,  0.3818300505051189, ...
          0.2797053914892766,  0.2797053914892766,  0.1294849661688697, ...
          0.1294849661688697];	
  case {14, 15} % R = 8
    Q = [-0.1834346424956498,  0.1834346424956498, -0.5255324099163290, ...
          0.5255324099163290, -0.7966664774136267,  0.7966664774136267, ...
         -0.9602898564975363,  0.9602898564975363]; 
    W = [ 0.3626837833783620,  0.3626837833783620,  0.3137066458778873, ...
          0.3137066458778873,  0.2223810344533745,  0.2223810344533745, ...
          0.1012285362903763,  0.1012285362903763];
  case {16, 17} % R = 9
    Q = [ 0.0000000000000000, -0.8360311073266358,  0.8360311073266358, ...
         -0.9681602395076261,  0.9681602395076261, -0.3242534234038089, ...
          0.3242534234038089, -0.6133714327005904,  0.6133714327005904];
    W = [ 0.3302393550012598,  0.1806481606948574,  0.1806481606948574, ...
          0.0812743883615744,  0.0812743883615744,  0.3123470770400029, ...
          0.3123470770400029,  0.2606106964029354,  0.2606106964029354];
  otherwise
    error('Order %d not implemented', qOrd)
end % switch
Q = (Q + 1)/2;  W = W/2; % transformation [-1; 1] -> [0, 1]
end % function
