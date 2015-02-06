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
function ret = assembleVecEdgePhiIntFuncCont(g, markE0Tbdr, funcCont, N)
global gPhi1D
p = (sqrt(8*N+1)-3)/2;
qOrd = 2*p+1; [Q, W] = quadRule1D(qOrd);
Q2X1   = @(X1,X2) g.B(:,1,1)*X1 + g.B(:,1,2)*X2 + g.coordV0T(:,1,1)*ones(size(X1));
Q2X2   = @(X1,X2) g.B(:,2,1)*X1 + g.B(:,2,2)*X2 + g.coordV0T(:,1,2)*ones(size(X1));
ret = zeros(g.numT, N);
for n = 1 : 3
  [Q1, Q2] = gammaMap(n, Q);
  cDn = funcCont(Q2X1(Q1, Q2), Q2X2(Q1, Q2));
  for i = 1 : N
    ret(:,i) = ret(:,i) + markE0Tbdr(:,n) .* ( cDn * (W' .* gPhi1D{qOrd}(:,i,n)) );
  end % for
end % for
ret = reshape(ret', g.numT*N, 1);
end % function
