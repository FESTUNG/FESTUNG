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
function [XP1, XP2] = theta(nn, np, X1, X2)
switch nn
  case 1
    switch np
      case 1,  XP1 = 1-X1;  XP2 = 1-X2;
      case 2,  XP1 = zeros(size(X1));  XP2 = X2;
      case 3,  XP1 = X1;  XP2 = zeros(size(X1));
    end % switch
  case 2
    switch np
      case 1,  XP1 = 1-X2;  XP2 = X2;
      case 2,  XP1 = zeros(size(X1));  XP2 = 1-X2;
      case 3,  XP1 = X2;  XP2 = zeros(size(X1));
    end % switch
  case 3
    switch np
      case 1,  XP1 = X1;  XP2 = 1-X1;
      case 2,  XP1 = zeros(size(X1));  XP2 = X1;
      case 3,  XP1 = 1-X1;  XP2 = zeros(size(X1));
    end % switch
end % switch
end % function
