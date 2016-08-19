% Evaluates a piecewise linear continuous function which is implicitly given by 
% its values in the vertices of a grid in a number of physical points on each 
% triangle.
%
%===============================================================================
%> @file evaluateFuncFromVertexValues.m
%>
%> @brief Evaluates a piecewise linear continuous function which is implicitly 
%>        given by its values in the vertices of a grid in a number of physical 
%>        points on each triangle.
%===============================================================================
%>
%> @brief Evaluates a piecewise linear continuous function which is implicitly 
%>        given by its values in the vertices of a grid in a number of physical 
%>        points on each triangle.
%>
%> This routine evaluates any piecewise linear and globally continuous function,
%> which is defined implicitly on a triangular grid by its values in each grid 
%> point.
%> Note that due to the backtransformation to the reference element and repeated
%> use of multiplications with bsxfun this method can be quite time consuming 
%> and is meant to be used during preprocessing only.
%>
%> @param g            The lists describing the geometric and topological 
%>                     properties of a triangulation (see 
%>                     <code>generateGridData()</code>)
%> @param vertexValues The values of the implicitly defined function in the 
%>                     vertices of the grid. @f$[numV \times 1]@f$
%>                     by <code>integrateRefEdgePhiIntPhiInt()</code>.
%>                     @f$[N \times N \times {N_\mathrm{data}}\times 3]@f$
%> @param X1           x1-coordinates of all evaluation points with each row
%>                     contaiuning the points located in the associated element
%>                     and the columns representing the numerous points per 
%>                     element.
%> @param X2           x2-coordinates of all evaluation points with each row
%>                     contaiuning the points located in the associated element
%>                     and the columns representing the numerous points per 
%>                     element.
%> @retval  cEval      A matrix containing all the evaluations of said function
%>                     with the eavluations arranged in the same way as X1 and
%>                     X2.
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
%>
function cEval = evaluateFuncFromVertexValues(g, vertexValues, X1, X2)

lenX = size(X1,2);
FinvX1 = 0.5 ./ (g.areaT * ones(1,lenX)) .* ( bsxfun(@times, g.B(:,2,2), X1) - bsxfun(@times, g.B(:,1,2), X2) - (g.B(:,2,2).*g.coordV0T(:,1,1) - g.B(:,1,2).*g.coordV0T(:,1,2) ) * ones(1,lenX) );
FinvX2 = 0.5 ./ (g.areaT * ones(1,lenX)) .* (-bsxfun(@times, g.B(:,2,1), X1) + bsxfun(@times, g.B(:,1,1), X2) + (g.B(:,2,1).*g.coordV0T(:,1,1) - g.B(:,1,1).*g.coordV0T(:,1,2) ) * ones(1,lenX) );
cEval  = bsxfun(@times, vertexValues(g.V0T(:,1)), 1 - FinvX1 - FinvX2) ...
       + bsxfun(@times, vertexValues(g.V0T(:,2)),     FinvX1         ) ...
       + bsxfun(@times, vertexValues(g.V0T(:,3)),             FinvX2 );
end % function
