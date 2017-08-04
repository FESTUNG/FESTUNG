% Evaluates the i-th basis function.

%===============================================================================
%> @file phi.m
%>
%> @brief Evaluates the i-th basis function.
%===============================================================================
%>
%> @brief Evaluates the @f$i@f$-th basis function on the reference triangle 
%>        @f$\hat{T}@f$ at points specified by a list of @f$\hat{x}^1@f$ 
%>        and @f$\hat{x}^2@f$ coordinates .
%>
%> @param  i   The index of the basis function.
%> @param  X1  A list of @f$\hat{x}^1@f$ coordinates.
%> @param  X2  A list of @f$\hat{x}^2@f$ coordinates.
%> @retval ret The @f$i@f$-th basis function in all points specified by
%>             <code>X1</code>, <code>X2</code>. It holds <code>size(X1) ==
%>             size(X2) == size(ret)</code>.
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
function ret = phi(i, X1, X2)
switch i
  case 1,  ret = sqrt(2)*ones(size(X1));
  case 2,  ret = 2 - 6*X1;
  case 3,  ret = 2*sqrt(3)*(1 - X1 - 2*X2);
  case 4,  ret = sqrt(6)*((10*X1 - 8).*X1 + 1);
  case 5,  ret = sqrt(3)*((5*X1 - 4).*X1 + (-15*X2 + 12).*X2 - 1);
  case 6,  ret = 3*sqrt(5)*((3*X1 + 8*X2 - 4).*X1 + (3*X2 - 4).*X2 + 1);
  case 7,  ret = 2*sqrt(2)*(-1+(15+(-45+35*X1).*X1).*X1); 
  case 8,  ret = 2*sqrt(6)*(-1+(13+(-33+21*X1).*X1).*X1+(2+(-24+42*X1).*X1).*X2);
  case 9,  ret = 2*sqrt(10)*(-1+(9+(-15+7*X1).*X1).*X1+(6+(-48+42*X1).*X1+(-6+42*X1).*X2).*X2);
  case 10, ret = 2*sqrt(14)*(-1+(3+(-3+X1).*X1).*X1+(12+(-24+12*X1).*X1+(-30+30*X1+20*X2).*X2).*X2);
  case 11, ret = sqrt(10)*(1+(-24+(126+(-224+126*X1).*X1).*X1).*X1);
  case 12, ret = sqrt(30)*(1+(-22+(105+(-168+84*X1).*X1).*X1).*X1+(-2+(42+(-168+168*X1).*X1).*X1).*X2);
  case 13, ret = 5*sqrt(2)*(1+(-18+(69+(-88+36*X1).*X1).*X1).*X1+(-6+(102+(-312+216*X1).*X1).*X1 ...
                   +(6+(-96+216*X1).*X1).*X2).*X2);
  case 14, ret = sqrt(70)*(1+(-12+(30+(-28+9*X1).*X1).*X1).*X1+(-12+(132+(-228+108*X1).*X1).*X1 ...
                   +(30+(-300+270*X1).*X1+(-20+180*X1).*X2).*X2).*X2);
  case 15, ret = 3*sqrt(10)*(1+(-4+(6+(-4+X1).*X1).*X1).*X1+(-20+(60+(-60+20*X1).*X1).*X1 ...
                   +(90+(-180+90*X1).*X1+(-140+140*X1+70*X2).*X2).*X2).*X2);
end % switch
end % function
