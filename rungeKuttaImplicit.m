% Provides weights and time levels for diagonally implicit Runge-Kutta methods 
% up to order 4.

%===============================================================================
%> @file rungeKuttaImplicit.m
%>
%> @brief Provides weights and time levels for diagonally implicit Runge-Kutta 
%> methods up to order 4.
%===============================================================================
%>
%> @brief Provides weights and time levels for diagonally implicit Runge-Kutta 
%> methods up to order 4.
%>
%> The present methods are all @f$A@f$- and @f$L@f$-stable. Furthermore, the 
%> methods are stiffly accurate, i.e. @f$b_j=a_{sj}@f$ where @f$s@f$ is the 
%> number of stages.
%> 
%> 
%> All other entries are zero.
%> @param  ord                TODO
%> 
%> @param  tau  TODO
%> 
%> @param  t0   TODO
%> 
%> 
%> @param  t    TODO
%> 
%> @param  A    TODO

%> @param  b    TODO
%> 
%> @retval c    TODO
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
function [t, A, b, c] = rungeKuttaImplicit(ord, tau, t0)
switch ord
  case 1
    A = 1;
    
  case 2
    lambda = 0.5 * (2 - sqrt(2));
    A = [     lambda,      0 ; 
          (1-lambda), lambda ];
    
  case 3
    alpha = 0.4358665215084589994160194511935568425292;
    beta  = 0.5 * (1.+ alpha);
    b1    = -0.25 * (6 * alpha * alpha - 16 * alpha + 1);
    b2    =  0.25 * (6 * alpha * alpha - 20 * alpha + 5);

    A = [        alpha,     0,     0 ; 
          beta - alpha, alpha,     0 ; 
                    b1,    b2, alpha ];
    
  case 4
    gamma = 1./4.; 
    A = [     gamma,           0,        0,        0,     0;
              1./2.,       gamma,        0,        0,     0;
            17./50.,     -1./25.,    gamma,        0,     0;
         371./1360., -137./2720., 15./544.,    gamma,     0;
            25./24.,    -49./48., 125./16., -85./12., gamma ];
    
  otherwise
    msgID = 'rungeKuttaImplicit:TimeIntegratorOrder';
    msg = 'Function is not defined for order<1 or order>4.';
    baseException = MException(msgID,msg);
    causeID = 'rungeKuttaImplicit:OrderTooLarge';
    causeMsg = sprintf('Value of order is %d', ord);
    causeException = MException(causeID,causeMsg);
    baseException = addCause(baseException,causeException);
    throw(baseException);
end % switch
    
b = A(end, :);
c = sum(A, 2);
t = t0 + c * tau;
end % function
