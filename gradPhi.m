% This file is part of FESTUNG 
% Copyright (C) 2014 Florian Frank, Balthasar Reuter, Vadym Aizinger
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
function ret = gradPhi(i, m, X1, X2)
switch m
  case 1
    switch i
      case 1,  ret = zeros(size(X1));
      case 2,  ret = -6*ones(size(X1));
      case 3,  ret = -2*sqrt(3)*ones(size(X1));
      case 4,  ret = sqrt(6)*(20*X1 - 8);
      case 5,  ret = sqrt(3)*(10*X1 - 4);
      case 6,  ret = 6*sqrt(5)*(3*X1 + 4*X2 - 2);
      case 7,  ret = 2*sqrt(2)*(15+(-90+105*X1).*X1);
      case 8,  ret = 2*sqrt(6)*(13+(-66+63*X1).*X1+(-24+84*X1).*X2);
      case 9,  ret = 2*sqrt(10)*(9+(-30+21*X1).*X1+(-48+84*X1+42*X2).*X2);
      case 10, ret = 2*sqrt(14)*(3+(-6+3*X1).*X1+(-24+24*X1+30*X2).*X2);
      case 11, ret = sqrt(10)*(-24+(252+(-672+504*X1).*X1).*X1);
      case 12, ret = sqrt(30)*(-22+(210+(-504+336*X1).*X1).*X1+(42+(-336+504*X1).*X1).*X2);
      case 13, ret = 5*sqrt(2)*(-18+(138+(-264+144*X1).*X1).*X1 ...
                       +(102+(-624+648*X1).*X1+(-96+432*X1).*X2).*X2);
      case 14, ret = sqrt(70)*(-12+(60+(-84+36*X1).*X1).*X1 ...
                       +(132+(-456+324*X1).*X1+(-300+540*X1+180*X2).*X2).*X2);
      case 15, ret = 3*sqrt(10)*(-4+(12+(-12+4*X1).*X1).*X1 ...
                       +(60+(-120+60*X1).*X1+(-180+180*X1+140*X2).*X2).*X2);
    end
  case 2
    switch i
      case 1,  ret = zeros(size(X1));
      case 2,  ret = zeros(size(X1));
      case 3,  ret = -4*sqrt(3)*ones(size(X1));
      case 4,  ret = zeros(size(X1));
      case 5,  ret = 2*sqrt(3)*( -15*X2 + 6);
      case 6,  ret = 6*sqrt(5)*(4*X1 + 3*X2 - 2);
      case 7,  ret = zeros(size(X1));
      case 8,  ret = 2*sqrt(6)*(2+(-24+42*X1).*X1);
      case 9,  ret = 2*sqrt(10)*(6+(-48+42*X1).*X1+(-12+84*X1).*X2);
      case 10, ret = 2*sqrt(14)*(12+(-24+12*X1).*X1+(-60+60*X1+60*X2).*X2);
      case 11, ret = zeros(size(X1));
      case 12, ret = sqrt(30)*(-2+(42+(-168+168*X1).*X1).*X1);
      case 13, ret = 5*sqrt(2)*(-6+(102+(-312+216*X1).*X1).*X1 ...
                       +(12+(-192+432*X1).*X1).*X2);
      case 14, ret = sqrt(70)*(-12+(132+(-228+108*X1).*X1).*X1 ...
                       +(60+(-600+540*X1).*X1+(-60+540*X1).*X2).*X2);
      case 15, ret = 3*sqrt(10)*(-20+(60+(-60+20*X1).*X1).*X1 ...
                       +(180+(-360+180*X1).*X1+(-420+420*X1+280*X2).*X2).*X2);
    end % switch
end % switch
end % function
