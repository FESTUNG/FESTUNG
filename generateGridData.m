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
function g = generateGridData(coordV, V0T)
g.coordV = coordV;
g.V0T = V0T;
g.numT = size(g.V0T, 1);
g.numV = size(g.coordV, 1);
% The following implicitely defines the signs of the edges.
g.V2T  = sparse(g.V0T(:, [1 2 3 1 2 3 1 2 3]), g.V0T(:, [2 3 1 2 3 1 2 3 1]), ...
  [(1:g.numT)',zeros(g.numT,3),(1:g.numT)',zeros(g.numT,3),(1:g.numT)'],g.numV,g.numV);
% The following implicitely defines the edge numbers.
[r, c] = find(triu(g.V2T + g.V2T'));
g.V2E = sparse(r, c, 1 : size(r, 1), g.numV, g.numV);
g.V2E = g.V2E + g.V2E';
idxE = full(g.V2E(sub2ind([g.numV,g.numV],g.V0T(end:-1:1,[1,2,3]),g.V0T(end:-1:1,[2,3,1]))))';
g.V0E(idxE(:), 1) = reshape(g.V0T(end:-1:1, [1,2,3])', 3*g.numT, 1);
g.V0E(idxE(:), 2) = reshape(g.V0T(end:-1:1, [2,3,1])', 3*g.numT, 1);
g.T0E(idxE(:), 1) = reshape(full(g.V2T(sub2ind([g.numV,g.numV], ...
  g.V0T(end:-1:1,[1,2,3]), g.V0T(end:-1:1,[2,3,1]))))', 3*g.numT, 1);
g.T0E(idxE(:), 2) = reshape(full(g.V2T(sub2ind([g.numV,g.numV], ...
  g.V0T(end:-1:1,[2,3,1]), g.V0T(end:-1:1,[1,2,3]))))', 3*g.numT, 1);
g.numE = size(g.V0E, 1);
vecE = g.coordV(g.V0E(:, 2), :) - g.coordV(g.V0E(:, 1), :);
g.areaE = (vecE(:, 1).^2 + vecE(:, 2).^2).^(1/2);
g.nuE = vecE * [0,-1; 1,0] ./ g.areaE(:, [1, 1]);
g.areaT = ...
  ( g.coordV(g.V0T(:,1),1).*g.coordV(g.V0T(:,2),2) + g.coordV(g.V0T(:,2),1).*g.coordV(g.V0T(:,3),2) ...
  + g.coordV(g.V0T(:,3),1).*g.coordV(g.V0T(:,1),2) - g.coordV(g.V0T(:,1),1).*g.coordV(g.V0T(:,3),2) ...
  - g.coordV(g.V0T(:,2),1).*g.coordV(g.V0T(:,1),2) - g.coordV(g.V0T(:,3),1).*g.coordV(g.V0T(:,2),2) )/2;
g.baryT = (g.coordV(g.V0T(:,1),:)+g.coordV(g.V0T(:,2),:)+g.coordV(g.V0T(:,3),:))/3;
g.E0T = full(g.V2E(sub2ind([g.numV,g.numV],g.V0T(:,[2,3,1]),g.V0T(:,[3,1,2]))));
g.areaE0T = g.areaE(g.E0T);
sigE0T = 1-2*(bsxfun(@eq, reshape(g.T0E(g.E0T,2),g.numT,3),(1:g.numT)'));
g.baryE = 0.5 * (g.coordV(g.V0E(:, 1), :) + g.coordV(g.V0E(:, 2), :));
for n = 1 : 3
  for m = 1 : 2
    g.coordV0T(:, n, m) = g.coordV(g.V0T(:, n), m)';
    g.baryE0T(:, n, m) = g.baryE(g.E0T(:, n), m)';
    g.nuE0T(:, n, m) = g.nuE(g.E0T(:, n), m)'.* sigE0T(:, n)';
  end 
  Tn = sigE0T(:, n) == 1;  Tp = ~Tn;
  g.E0E(g.E0T(Tn, n), 1) = n;  g.E0E(g.E0T(Tp, n), 2) = n;
end % for
for m = 1 : 2
  g.B(:, m, 1) = g.coordV0T(:, 2, m) - g.coordV0T(:, 1, m);
  g.B(:, m, 2) = g.coordV0T(:, 3, m) - g.coordV0T(:, 1, m);
end % for
markEint = g.E0E(:, 2) ~= 0; % mark interior edges
g.markE0TE0T = cell(3, 3);
for nn = 1 : 3
  for np = 1 : 3
    g.markE0TE0T{nn,np} = sparse(g.numT, g.numT);
    markEn = g.E0E(:, 1) == nn;  markEp = g.E0E(:, 2) == np;
    idx = markEn & markEp & markEint;
    g.markE0TE0T{nn, np}(sub2ind([g.numT, g.numT], g.T0E(idx, 1), g.T0E(idx, 2))) = 1;
    markEn = g.E0E(:, 2) == nn;  markEp = g.E0E(:, 1) == np;
    idx = markEn & markEp & markEint;
    g.markE0TE0T{nn, np}(sub2ind([g.numT, g.numT], g.T0E(idx, 2), g.T0E(idx, 1))) = 1;
  end % for
end % for
end % function
