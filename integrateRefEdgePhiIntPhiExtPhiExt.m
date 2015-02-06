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
function ret = integrateRefEdgePhiIntPhiExtPhiExt(N)
global gPhi1D gThetaPhi1D
p = (sqrt(8*N+1)-3)/2;  qOrd = 2*p+1;  [~, W] = quadRule1D(qOrd);
ret = zeros(N,N,N,3,3); % [N x N x N x 3 x 3]
for nn = 1 : 3 % 3 edges
  for np = 1 : 3
    for l = 1 : N
      for i = 1 : N
        for j = 1 : N
          ret(i, j, l, nn,np) = sum( W'.*gPhi1D{qOrd}(:,i,nn) .* ...
            gThetaPhi1D{qOrd}(:,l,nn,np) .* gThetaPhi1D{qOrd}(:,j,nn,np) );
        end % for
      end % for
    end % for
  end % for
end % for
end
