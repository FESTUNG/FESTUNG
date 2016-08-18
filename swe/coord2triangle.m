% Uses the evaluation of the unknowns of the 2D Shallow-Water Equations to 
% compute their harmonic averages of Roe-Pike type for interior edges.
%
%===============================================================================
%> @file computeAveragedVariablesQ0E0Tint.m
%>
%> @brief Uses the evaluation of the unknowns of the 2D Shallow-Water Equations 
%>        to compute their harmonic averages of Roe-Pike type for interior edges
%>        (see @ref ROE1985 for details)
%===============================================================================
%>
%> @brief Uses the evaluation of the unknowns of the 2D Shallow-Water Equations 
%>        to compute their harmonic averages of Roe-Pike type for interior edges
%>        (see @ref ROE1985 for details)
%>
%> This routine uses the limits of the height and momenta on both sides of the
%> edges to compute a stable avarage, later to be used for flux approximation.
%> There are different types of averages available:
%> full-harmonic uses the Roe-Pike approach on both the normal component of the
%> velocity and the root of the graviational constant times the height.
%> semi-harmonic uses this approach on the normal component of the
%> velocity and takes the mean of the height for the rest.
%> mean uses mean values everywhere.
%>
%> @param g           The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%> @param coord1      List of physical x1-coordinates of stations.
%> @param coord2      List of physical x1-coordinates of stations
%>
%> @retval ret        Cell that for each element stores the triangles and the 
%>                    barycentric coordinates of the corresponding station.
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Balthasar Reuter, Florian Frank, Vadym Aizinger
%>                      Modified 08/18/16 by Hennes Hajduk
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
function ret = coord2triangle(g, coord1, coord2)
numCoord = length(coord1);
% Compute barycentric coordinates for each cartesian coordinate pair
coordBary = zeros(numCoord, g.numT, 3);
coordBary(:,:,1) = bsxfun(@times, (g.coordV0T(:,2,2) - g.coordV0T(:,3,2)).', ...
                          bsxfun(@minus, coord1, g.coordV0T(:,3,1).')) + ...
                   bsxfun(@times, (g.coordV0T(:,3,1) - g.coordV0T(:,2,1)).', ...
                          bsxfun(@minus, coord2, g.coordV0T(:,3,2).'));
coordBary(:,:,2) = bsxfun(@times, (g.coordV0T(:,3,2) - g.coordV0T(:,1,2)).', ...
                          bsxfun(@minus, coord1, g.coordV0T(:,3,1).')) + ...
                   bsxfun(@times, (g.coordV0T(:,1,1) - g.coordV0T(:,3,1)).', ...
                          bsxfun(@minus, coord2, g.coordV0T(:,3,2).'));
coordBary(:,:,1:2) = bsxfun(@times, coordBary(:,:,1:2), .5 ./ g.areaT.');
coordBary(:,:,3) = 1 - coordBary(:,:,1) - coordBary(:,:,2);
% Find triangles with only positive barycentric coordinates
[idxCoord, idxTriangle] = find(coordBary(:,:,1) >= 0 & coordBary(:,:,2) >= 0 & coordBary(:,:,3) >= 0);
% TODO: vectorize this
ret = cell(numCoord,1);
for n = 1 : numCoord
  triangles = idxTriangle(idxCoord == n);
  ret{n} = [triangles, reshape(coordBary(n,triangles,:), length(triangles), 3) ];
end % for
end % function
