% Routine that uses edge bisection of triangles to refine any triangular grid as
% given by any of the grid generating routines in FESTUNG.

%===============================================================================
%> @file refineGrid.m
%>
%> @brief Routine that uses edge bisection of triangles to refine any triangular  
%>        grid as as given by any of the grid generating routines in FESTUNG.
%===============================================================================
%>
%> @brief Routine that uses edge bisection of triangles to refine any triangular 
%>        grid as given by any of the grid generating routines in FESTUNG.
%>
%> @param  g            A struct with all numerical and logical fields or cells
%>                      identifying the mesh of a compuational domain that are
%>                      required in FESTUNG as computed by 
%>                      <code>generateGridData()</code> and identifying boundary
%>                      edge fields idE and idE0T.
%>
%> @retval g            The struct that is the result of refining the input 
%>                      grid by edge bisection of each triangle faeturing all
%>                      inforamtion provided by <code>generateGridData()</code>
%>                      and idE and idE0T.
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Hennes Hajduk, Balthasar Reuter, Florian Frank, 
%>                      Vadym Aizinger
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
function [g, dataV] = refineGrid(g, dataL)
K = g.numT;
vertCtr = g.numV;
edges = zeros(g.numE,1);
V0T = zeros(4*K,3);
coordV = [g.coordV; zeros(g.numE,2)];
idE0T = zeros(4*K,3);
bdrTypes = max(g.idE);
dataV = zeros(g.numV,1);
for k = 1 : K
	% old vertices
	V0T(4*(k-1)+2,1) = g.V0T(k,1);
	V0T(4*(k-1)+3,2) = g.V0T(k,2);
	V0T(4*(k-1)+4,3) = g.V0T(k,3);
  dataV(g.V0T(k,1)) = dataL(k,1);
  dataV(g.V0T(k,2)) = dataL(k,2);
  dataV(g.V0T(k,3)) = dataL(k,3);
	% local edge 1
	edge = g.E0T(k,1);
	if edges(edge) ~= 0
		V0T(4*(k-1)+1,1) = edges(edge);
		V0T(4*(k-1)+3,3) = edges(edge);
		V0T(4*(k-1)+4,2) = edges(edge);
	else
		vertCtr	= vertCtr+1;
		edges(edge)	= vertCtr;
		coordV(vertCtr,:) = 0.5*(coordV(g.V0E(edge,1),:)+coordV(g.V0E(edge,2),:));
		V0T(4*(k-1)+1,1) = vertCtr;
		V0T(4*(k-1)+3,3) = vertCtr;
		V0T(4*(k-1)+4,2) = vertCtr;
    
    dataV(vertCtr) = 0.5*(dataL(k,2) + dataL(k,3));
    
		% boundary edges
		if g.idE(edge) ~= 0
			idE0T(4*(k-1)+3,1) = g.idE(edge);
			idE0T(4*(k-1)+4,1) = g.idE(edge);
		end % if
	end % if
	% local edge 2
	edge = g.E0T(k,2);
	if edges(edge) ~= 0
		V0T(4*(k-1)+1 ,2) = edges(edge);
		V0T(4*(k-1)+2 ,3) = edges(edge);
		V0T(4*(k-1)+4 ,1) = edges(edge);
	else
		vertCtr	= vertCtr+1;
		edges(edge)	= vertCtr;
		coordV(vertCtr,:)	= 0.5*(coordV(g.V0E(edge,1),:)+coordV(g.V0E(edge,2),:));
		V0T(4*(k-1)+1,2) = vertCtr;
		V0T(4*(k-1)+2,3) = vertCtr;
		V0T(4*(k-1)+4,1) = vertCtr;
    
    dataV(vertCtr) = 0.5*(dataL(k,3) + dataL(k,1));
    
		% boundary edges
		if g.idE(edge) ~= 0
			idE0T(4*(k-1)+2,2) = g.idE(edge);
			idE0T(4*(k-1)+4,2) = g.idE(edge);
		end % if
	end % if
	% local edge 3
	edge = g.E0T(k,3);
	if edges(edge) ~= 0
		V0T(4*(k-1)+1,3) = edges(edge);
		V0T(4*(k-1)+2,2) = edges(edge);
		V0T(4*(k-1)+3,1) = edges(edge);
	else
		vertCtr	= vertCtr+1;
		edges(edge)	= vertCtr;
		coordV(vertCtr,:)	= 0.5*(coordV(g.V0E(edge,1),:)+coordV(g.V0E(edge,2),:));
		V0T(4*(k-1)+1,3) = vertCtr;
		V0T(4*(k-1)+2,2) = vertCtr;
		V0T(4*(k-1)+3,1) = vertCtr;
    
    dataV(vertCtr) = 0.5*(dataL(k,1) + dataL(k,2));
    
		% boundary edges
		if g.idE(edge) ~= 0
			idE0T(4*(k-1)+2,3) = g.idE(edge);
			idE0T(4*(k-1)+3,3) = g.idE(edge);
		end % if
	end % if
end % for
g = generateGridData(coordV,V0T);
g.idE = zeros(g.numE,1);
for b = 1:bdrTypes
  g.idE(g.E0T(idE0T==b)) = b;
end % for
g.idE0T = idE0T;
end % function