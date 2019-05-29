% Rotates a grid struct around zero by a certain degree.
%
%===============================================================================
%> @file rotateGrid.m
%>
%> @brief Rotates a grid struct around zero by a certain degree.
%===============================================================================
%>
%> @brief Rotates a grid struct around zero by a certain degree.
%>
%> The specified angle must be in circular measure. 
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of the unrotated triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  angle      The rotation angle to be used.
%> @retval g          The lists describing the geometric and topological 
%>                    properties of the rotated triangulation
%>                    @f$[1 \times 1 \text{ struct}]@f$
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
function g = rotateGrid(g, angle)

idE = g.idE;

c = cos(angle);
s = sin(angle);

coordV = [c * g.coordV(:,1)-s * g.coordV(:,2), s * g.coordV(:,1)+c * g.coordV(:,2)];

g = generateGridData(coordV, g.V0T);

g.idE = idE;
g.idE0T = g.idE(g.E0T);

end % function
