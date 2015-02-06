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
function dataLagr = projectDataDisc2DataLagr(dataDisc)
[K, N] = size(dataDisc);
switch N
  case 1,    L1 = 1/3;                     L2 = 1/3;                    % locally constant
  case 3,    L1 = [0, 1, 0];               L2 = [0, 0, 1];              % locally linear
  otherwise, L1 = [0, 1, 0, 1/2, 0, 1/2];  L2 = [0, 0, 1, 1/2, 1/2, 0]; % locally quadratic
end % switch
dataLagr = zeros(K, length(L1));
for i = 1 : N
  dataLagr = dataLagr + dataDisc(:, i) * phi(i, L1, L2);
end % for
end % function
