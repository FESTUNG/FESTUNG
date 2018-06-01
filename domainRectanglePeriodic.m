% Creates a structured grid grid over a rectangular domain with periodic 
% boundary conditions.
%
%===============================================================================
%> @file
%>
%> @brief Creates a structured grid grid over a rectangular domain with periodic
%>        boundary conditions. South and north boundaries as well as east
%>        and west boundaries can be associated with each other.
%===============================================================================
%>
%> @brief Creates a structured grid grid over a rectangular domain with periodic
%>        boundary conditions. South and north boundaries as well as east
%>        and west boundaries can be associated with each other.
%>
%> First a structured grid is created. For the purpose of bandwidth it makes 
%> sense to have fewer elements in y- than in x-direction. 
%> Then the corresponding edges and corners are used to modify the
%> topological grid data using the routine <code>modifyGridDataPeriodic()<code>
%> in order to prescribe the desired periodic boundary conditions.
%>
%> @param  g              The lists describing the geometric and topological 
%>                        properties of a triangulation (see 
%>                        <code>generateGridData()</code>) 
%>                        @f$[1 \times 1 \text{ struct}]@f$
%> @param  X              The x-coordinates of the rectangle corners
%> @param  Y              The y-coordinates of the rectangle corners
%> @param  hx             The length of the x-cathetus of each element
%> @param  hy             The length of the y-cathetus of each element
%> @param  identification A string indicating which boundaries to identify.
%>                        WE-SN identifies the left and right as well as
%>                        lower and upper boundaries. WE identifies the 
%>                        left and right boundaries. SN identifies the
%>                        lower and upper boundaries.
%> @retval g              The modified lists describing the geometric and
%>                        topological properties of a triangulation (see 
%>                        <code>generateGridData()</code>) 
%>                        @f$[1 \times 1 \text{ struct}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2018 Hennes Hajduk, Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function g = domainRectanglePeriodic(X, Y, hx, hy, identification)

if nargin < 4
  hy = hx;
  identification = 'WE-SN';
elseif nargin < 5
  identification = 'WE-SN';
end % if

[g, isPeriodicMarkV0TV0T] = domainRectangle(X, Y, hx, hy, true, identification);
numTy = length(Y(1):hy:Y(2)) - 1;

if strcmp(identification, 'WE-SN')
  periodicEdges = [ g.E0T(1:2:2*numTy-1,2) g.E0T(end-2*numTy+2:2:end,3); ...  % west and east
                    g.E0T(1:2*numTy:end,3) g.E0T(2*numTy:2*numTy:end,1) ];    % south and north

  cornerInd = [1, numTy+1, g.numV-numTy g.numV];
elseif strcmp(identification, 'WE')
  periodicEdges = [ g.E0T(1:2:2*numTy-1,2) g.E0T(end-2*numTy+2:2:end,3) ];    % west and east
  cornerInd = [];
elseif strcmp(identification, 'SN')
  periodicEdges = [ g.E0T(1:2*numTy:end,3) g.E0T(2*numTy:2*numTy:end,1) ];    % south and north
  cornerInd = [];
else
  error('Unsupported combination of periodic boundary edges.')
end % if

g = modifyGridDataPeriodic(g, periodicEdges, cornerInd, isPeriodicMarkV0TV0T);
end % function
