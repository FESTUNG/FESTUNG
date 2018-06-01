% Generates a union-jack type mesh for a square.

%===============================================================================
%> @file
%>
%> @brief Generates a union-jack type mesh for a square.
%===============================================================================
%>
%> @brief Generates a union-jack type mesh for a square.
%>
%> The four corners, given by (X(1), Y(1)), (X(2), Y(1)), (X(2), Y(2)), (X(1), Y(2)) 
%> form a rectangle, that is subdivided into triangles according to the mesh 
%> parameters hx and hy for each spatial dimension. The subroutine stencil
%> is called to compute the indices of the element vertices according to whether
%> the current center of the union jack is next to an edge or not.
%> In order to simplify the indexing a counter variable is handed back and forth
%> between the routines domainUnionJack and stencil.
%>
%> Note that this routine is not vectorized.
%>
%> @param  X     The two x-coordinates of the corners of the domain.
%> @param  Y     The two y-coordinates of the corners of the domain.
%> @param  hx   The mesh parameter in x-direction.
%> @param  hy   The mesh parameter in y-direction.
%> @retval g       A struct containing the lists that describe the
%>                triangulation, as explained in <code>generateGridData()</code>.
%>                <code>idE</code> and <code>ideE0T</code> are filled correctly.
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2018 Florian Frank, Balthasar Reuter, Vadym Aizinger
%> 
%> @author Hennes Hajduk, 2018.
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
function g = domainUnionJack(X, Y, hx, hy)

if nargin < 4
  hy = hx;
end % if

x = X(1):hx:X(2);
y = Y(1):hy:Y(2);

numTx = length(x)-1;
numTy = length(y)-1;

[X, Y] = meshgrid(x,y);
coordV = [reshape(X, [], 1), reshape(Y, [], 1)];
V0T = zeros(2*numTx*numTy, 3);

counter = 1;
for i = 2:2:numTx+1
	for j = 2:2:numTy+1
    [V0T, counter] = stencil(V0T, i, j, numTx, numTy, counter);
	end
end
g = generateGridData(coordV, V0T);
g.idE = zeros(g.numE,1);

g.idE(g.baryE(:, 2) == Y(1)) = 1; % south
g.idE(g.baryE(:, 1) == X(2)) = 2; % east
g.idE(g.baryE(:, 2) == Y(2)) = 3; % north
g.idE(g.baryE(:, 1) == X(1)) = 4; % west

g.idE0T = g.idE(g.E0T);
end % function 

function [V0T, counter] = stencil(V0T, i, j, numTx, numTy, counter)
centerInd = (i-1)*(numTy+1)+j;
if (i < numTx+1)
  if (j < numTy+1)
    indList = [ centerInd-numTy-2 centerInd-1 centerInd-numTy-1;
                centerInd-numTy-1 centerInd+1 centerInd-numTy;
                centerInd-numTy-1 centerInd-1 centerInd;
                centerInd-numTy-1 centerInd centerInd+1;
                centerInd centerInd-1 centerInd+numTy+1;
                centerInd centerInd+numTy+1 centerInd+1;
                centerInd-1 centerInd+numTy centerInd+numTy+1;
                centerInd+1 centerInd+numTy+1 centerInd+numTy+2 ];
    V0T(counter:counter+7,:) = indList;
    counter = counter + 8;
  else
    indList = [ centerInd-numTy-2 centerInd-1 centerInd-numTy-1;
                centerInd-numTy-1 centerInd-1 centerInd;
                centerInd centerInd-1 centerInd+numTy+1;
                centerInd-1 centerInd+numTy centerInd+numTy+1 ];
    V0T(counter:counter+3,:) = indList;
    counter = counter + 4;
  end % if
else
  if (j < numTy+1)
    indList = [ centerInd-numTy-2 centerInd-1 centerInd-numTy-1;
                centerInd-numTy-1 centerInd+1 centerInd-numTy;
                centerInd-numTy-1 centerInd-1 centerInd;
                centerInd-numTy-1 centerInd centerInd+1 ];
    V0T(counter:counter+3,:) = indList;
    counter = counter + 4;
  else
    indList = [ centerInd-numTy-2 centerInd-1 centerInd-numTy-1;
                centerInd-numTy-1 centerInd-1 centerInd; ];
    V0T(counter:counter+1,:) = indList;
    counter = counter + 2;
  end % if
end % if
end % function
