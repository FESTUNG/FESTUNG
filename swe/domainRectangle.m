% Creates a structured grid grid over a rectangular domain.
%
%===============================================================================
%> @file
%>
%> @brief Creates a structured grid grid over a rectangular domain.
%===============================================================================
%>
%> @brief Creates a structured grid grid over a rectangular domain.
%>
%> A structured grid is created. For the purpose of bandwidth it makes 
%> sense to have fewer elements in y- than in x-direction. 
%>
%> @param  g                     The lists describing the geometric and topological 
%>                               properties of a triangulation (see 
%>                               <code>generateGridData()</code>) 
%>                               @f$[1 \times 1 \text{ struct}]@f$
%> @param  X                     The x-coordinates of the rectangle corners
%> @param  Y                     The y-coordinates of the rectangle corners
%> @param  hx                    The length of the x-cathetus of each element
%> @param  hy                    The length of the y-cathetus of each element
%> @param  isPeriodic            Logical that is used to decide whether to impose
%>                               periodic boundary conditions.
%> @param  identification        A string indicating which boundaries to identify.
%>                               WE-SN identifies the left and right as well as
%>                               lower and upper boundaries. WE identifies the 
%>                               left and right boundaries. SN identifies the
%>                               lower and upper boundaries.
%> @retval g                     The modified lists describing the geometric and
%>                               topological properties of a triangulation (see 
%>                               <code>generateGridData()</code>) 
%>                               @f$[1 \times 1 \text{ struct}]@f$
%> @retval isPeriodicMarkV0TV0T  Logical that indicates if the field
%>                               markV0TV0T is computed using a slow but
%>                               memory efficient algorithm to account for
%>                               periodic boundaries.
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2018 Florian Frank, Balthasar Reuter, Vadym Aizinger
%>
%> @author Hennes Hajduk, 2018
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
function [g, isPeriodicMarkV0TV0T] = domainRectangle(X, Y, hx, hy, isPeriodic, identification)

if nargin < 4
  hy = hx;
  isPeriodic = false;
  identification = 'WE-SN';
elseif nargin < 5
  isPeriodic = false;
  identification = 'WE-SN';
elseif nargin < 6
  identification = 'WE-SN';
end % if

assert(isPeriodic || strcmp(identification, 'WE-SN'), 'Boundary edge identifiaction is only possible for periodic boundaries.')

x = X(1):hx:X(2);
y = Y(1):hy:Y(2);

numTx = length(x)-1;
numTy = length(y)-1;

[X, Y] = meshgrid(x,y);
coordV = [reshape(X, [], 1), reshape(Y, [], 1)];
V0T1 = zeros(3, numTx*numTy);
V0T2 = zeros(3, numTx*numTy);

for i = 1:numTx
	for j = 1:numTy
		V0T1(:,(i-1)*numTy+j) = [j j+numTy+1 j+1] + (i-1)*(numTy+1);
		V0T2(:,(i-1)*numTy+j) = [j+numTy+1 j+numTy+2 j+1] + (i-1)*(numTy+1);
	end
end
V0T = reshape([V0T1; V0T2] , 3, [])';
g = generateGridData(coordV, V0T);

if ~isfield(g, 'markV0TV0T')
  if ~strcmp(identification, 'WE-SN')
    error('So far it is only possible to have periodic boundaries on all boundary edges.') % TODO
  end % if
  warning('Usual computation of sparse field markV0TV0T required too much memory. Using memory efficient but non vectorized version instead. Long preprocessing time to be expected.')
  g.markV0TV0T = computeMarkV0TV0T(numTx, numTy, isPeriodic);
  isPeriodicMarkV0TV0T = true;
else
  isPeriodicMarkV0TV0T = false;
end % if

g.idE = zeros(g.numE,1);

g.idE(g.baryE(:, 2) == y(1))   = 1; % south
g.idE(g.baryE(:, 1) == x(end)) = 2; % east
g.idE(g.baryE(:, 2) == y(end)) = 3; % north
g.idE(g.baryE(:, 1) == x(1))   = 4; % west

g.idE0T = g.idE(g.E0T);

end % function
