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
function g = domainPolygon(X1, X2, h)
gd = [2; length(X1(:)); X1(:); X2(:)]; % geometry description
sf = 'polygon';                        % set formula
ns = double('polygon')';               % name space
[p, e, t] = initmesh(decsg(gd,sf,ns), 'Hmax', h);
g  = generateGridData(p', t(1:3, :)');
g.idE = zeros(g.numE, 1);
g.idE(g.V2E(sub2ind(size(g.V2E),e(1,:),e(2,:)))) = e(5,:);
g.idE0T = g.idE(g.E0T); % local edge IDs
end % function
