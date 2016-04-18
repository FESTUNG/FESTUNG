% Generates a mesh for the unit circle with given mesh width.
%
%===============================================================================
%> @file domainCircle.m
%>
%> @brief Generates a mesh for the unit circle with given mesh width.
%===============================================================================
%>
%> @brief Generates a mesh for the unit circle with given mesh width.
%>
%> This function calls Gmsh with the geometry defined in 
%> <code>domainCircle.geo</code>:
%> @verbinclude domainCircle.geo
%> and the maximum element width <code>h</code>.
%> 
%> Gmsh produces a file <code>domainCircle.mesh</code>, which is used to 
%> generate all required lists and call <code>generateGridData()</code>.
%>
%> @par Example
%> @parblock
%> @code
%> g = domainCircle(1/3);
%> visualizeGrid(g)
%> @endcode
%> produces the following output:
%> @image html  domainCircle.png
%> @endparblock
%>
%> @param  h  The maximum diameter of an element.
%> @retval g  A struct containing the lists that describe the
%>            triangulation, as explained in <code>generateGridData()</code>.
%>            <code>idE</code> and <code>ideE0T</code> are filled correctly.
%>
%> @note  Gmsh (http://geuz.org/gmsh/) must be installed and available in the
%>        <code>PATH</code>-variable. If MATLAB can't find Gmsh, add Gmshs path  
%>        to the variable:
%>        @code
%>        setenv('PATH', [getenv('PATH'), ':/path/to/gmsh/bin'])
%>        @endcode
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2015 Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function g = domainCircle(h)
%% Generation of domainCircle.mesh using domainCircle.geo.
cmd = sprintf('gmsh -2 -format mesh -clscale %f -o "domainCircle.mesh" "domainCircle.geo"' , h);
system(cmd);
%% Extract data from domainCircle.mesh.
fid = fopen('domainCircle.mesh', 'r');
tline = fgets(fid);
while ischar(tline)
  if strfind(tline, 'Vertices')
    numV = fscanf(fid, '%d', [1, 1]);
    coordV = reshape(fscanf(fid, '%f'), 4, numV)';  coordV(:, 3:4) = [];
  end % if
  if strfind(tline, 'Edges')
    numEbdry = fscanf(fid, '%d', [1, 1]);
    tmp = reshape(fscanf(fid, '%f'), 3, numEbdry)';
    V0Ebdry = tmp(:, 1:2);  idEbdry = tmp(:, 3);
  end % if
  if strfind(tline, 'Triangles')
    numT = fscanf(fid, '%d', [1, 1]);
    V0T = reshape(fscanf(fid, '%f'), 4, numT)';  V0T(:, 4) = [];
  end % if
  tline = fgets(fid);
end % while
fclose(fid);
%% Generate lists and set boundary IDs.
g = generateGridData(coordV, V0T);
g.idE = zeros(g.numE, 1);
g.idE(g.V2E(sub2ind([g.numV,g.numV], V0Ebdry(:,1), V0Ebdry(:,2)))) = idEbdry;
g.idE0T = g.idE(g.E0T); % local edge IDs
end % function
