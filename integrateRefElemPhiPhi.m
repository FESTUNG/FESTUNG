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
function ret = integrateRefElemPhiPhi(N)
global gPhi2D
p = (sqrt(8*N+1)-3)/2;  qOrd = max(2*p, 1);  [~,~,W] = quadRule2D(qOrd);
ret = zeros(N); % [N x N]
for i = 1 : N
  for j = 1 : N
    ret(i, j) = sum( W' .* gPhi2D{qOrd}(:, i) .* gPhi2D{qOrd}(:, j) );
  end % for
end % for
end % function
