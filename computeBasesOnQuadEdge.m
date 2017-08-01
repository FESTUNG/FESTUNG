% TODO
%===============================================================================
%> @file computeBasesOnQuadEdge.m
%>
%> @brief TODO
%===============================================================================
%>
%> @brief TODO
%>
%> TODO
%> 
%> All other entries are zero.
%> @param  N                TODO
%> 
%> @param  basesOnQuadEdge  TODO
%> 
%> @param  requiredOrders   TODO
%> 
%> 
%> @retval ret              TODO
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
%> @author Balthasar Reuter, 2017
%> 
%> @par License
%> @parblock
%> This program is free software: you can redistribute it and/or modify
%> it under the terms of the GNU General Public License as published by
%> the Free Software Foundation, either version 3 of the License, or
%> (at your option) any later version.
%>
%> This program is distributed in the hope that it will be useful,
%> but WITHOUT ANY WARRANTY; without even the implied warranty of
%> MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%> GNU General Public License for more details.
%>
%> You should have received a copy of the GNU General Public License
%> along with this program.  If not, see <http://www.gnu.org/licenses/>.
%> @endparblock
%
function basesOnQuadEdge = computeBasesOnQuadEdge(N, basesOnQuadEdge, requiredOrders)
validateattributes(basesOnQuadEdge, {'struct'}, {}, mfilename, 'basesOnQuadEdge')

if nargin < 3
    p = N - 1;
    if p > 0
        requiredOrders = [2*p, 2*p+1];
    else
        requiredOrders = 1;
    end % if
end % if

% Precompute basis functions
basesOnQuadEdge.mu = cell(max(requiredOrders),1);
basesOnQuadEdge.thetaMu = cell(max(requiredOrders),1);
for it = 1 : length(requiredOrders)
    ord = requiredOrders(it);
    [Q, ~] = quadRule1D(ord);
    R = length(Q);
    
    basesOnQuadEdge.mu{ord} = zeros(R, N);
    for i = 1 : N
        basesOnQuadEdge.mu{ord}(:, i) = phi1D(i, Q);
    end
    
    basesOnQuadEdge.thetaMu{ord} = zeros(R, N, 2);
    basesOnQuadEdge.thetaMu{ord}(:, :, 1) = basesOnQuadEdge.mu{ord};
    basesOnQuadEdge.thetaMu{ord}(:, :, 2) = flipud(basesOnQuadEdge.mu{ord});
end % for
end % function
