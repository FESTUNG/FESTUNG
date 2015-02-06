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
function ret = assembleMatEdgePhiPhiFuncDiscNu(g, markE0Tint, refEdgePhiIntPhiIntPhiInt, refEdgePhiIntPhiExtPhiExt, dataDisc)
[K, N] = size(dataDisc);  
ret = cell(2, 1); ret{1} = sparse(K*N, K*N); ret{2} = sparse(K*N, K*N);
for nn = 1 : 3
  Rkn = 0.5 * g.areaE0T(:,nn);
  for np = 1 : 3
    markE0TE0TtimesRkn1 = bsxfun(@times, g.markE0TE0T{nn, np}, Rkn.*g.nuE0T(:,nn,1));
    markE0TE0TtimesRkn2 = bsxfun(@times, g.markE0TE0T{nn, np}, Rkn.*g.nuE0T(:,nn,2));
    for l = 1 : N
      ret{1} = ret{1} + kron(bsxfun(@times, markE0TE0TtimesRkn1, dataDisc(:,l).'), refEdgePhiIntPhiExtPhiExt(:,:,l,nn,np));
      ret{2} = ret{2} + kron(bsxfun(@times, markE0TE0TtimesRkn2, dataDisc(:,l).'), refEdgePhiIntPhiExtPhiExt(:,:,l,nn,np));
    end % for
  end % for
  Rkn = Rkn .* markE0Tint(:, nn);
  for l = 1 : N
    ret{1} = ret{1} + kron(spdiags(Rkn.*g.nuE0T(:,nn,1).*dataDisc(:,l),0,K,K), refEdgePhiIntPhiIntPhiInt(:,:,l,nn));
    ret{2} = ret{2} + kron(spdiags(Rkn.*g.nuE0T(:,nn,1).*dataDisc(:,l),0,K,K), refEdgePhiIntPhiIntPhiInt(:,:,l,nn));
  end % for
end % for
end % function