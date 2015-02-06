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
function err = computeL2Error(g, dataDisc, funcCont, qOrd)
global gPhi2D
N = size(dataDisc, 2); qOrd = max(qOrd,1);
[Q1, Q2, W] = quadRule2D(qOrd);
R = length(W);
X1 = kron(g.B(:,1,1),Q1)+kron(g.B(:,1,2),Q2)+kron(g.coordV0T(:,1,1),ones(1,R));
X2 = kron(g.B(:,2,1),Q1)+kron(g.B(:,2,2),Q2)+kron(g.coordV0T(:,1,2),ones(1,R));
cExOnQuadPts = funcCont(X1, X2); % [K x R]
cApprxOnQuadPts = dataDisc*gPhi2D{qOrd}'; % [K x R] = [K x N] * [N x R]
err = sqrt(2*dot((cApprxOnQuadPts - cExOnQuadPts).^2 * W.', g.areaT)); 
end % function
