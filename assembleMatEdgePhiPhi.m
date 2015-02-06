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
function ret = assembleMatEdgePhiPhi(g, markE0Tint, refEdgePhiIntPhiInt, refEdgePhiIntPhiExt)
K = g.numT;  N = size(refEdgePhiIntPhiInt, 1);
ret = sparse(K*N, K*N);
for n = 1 : 3
  ret = ret + kron(spdiags(markE0Tint(:,n),0,K,K), refEdgePhiIntPhiInt(:,:,n));
end % for
for nn = 1 : 3
  for np = 1 : 3
    ret = ret - kron(g.markE0TE0T{nn, np}, refEdgePhiIntPhiExt(:,:,nn,np));
  end % for
end % for
end % function
