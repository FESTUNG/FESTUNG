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
function dataDisc = projectFuncCont2DataDisc(g, funcCont, ord, refElemPhiPhi)
global gPhi2D
ord = max(ord,1);  [Q1, Q2, W] = quadRule2D(ord);
K = g.numT; N = size(refElemPhiPhi, 1);
F1 = @(X1, X2) g.B(:,1,1)*X1 + g.B(:,1,2)*X2 + g.coordV0T(:,1,1)*ones(size(X1));
F2 = @(X1, X2) g.B(:,2,1)*X1 + g.B(:,2,2)*X2 + g.coordV0T(:,1,2)*ones(size(X1));
rhs = funcCont(F1(Q1, Q2), F2(Q1, Q2)) * (repmat(W.', 1, N) .* gPhi2D{ord});
dataDisc = rhs / refElemPhiPhi;
end % function
