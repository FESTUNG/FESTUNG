% Reads grid and boundary parameters from fort.14- and fort.17-files and creates
% a grid as well as fields that store model parameter data.

%===============================================================================
%> @file
%>
%> @brief Reads grid and boundary parameters from fort.14- and fort.17-files and
%>        creates a grid as well as fields that store model parameter data.
%===============================================================================
%>
%> @brief Reads grid and boundary parameters from fort.14- and fort.17-files and
%>        creates a grid as well as fields that store model parameter data.
%>
%> This function reads in the coordinates, connectivity and values for the depth
%> in each vertex of a grid from the fort.14-file, see
%> <http://www.unc.edu/ims/adcirc/documentv47/fort_14.html>. This is used to 
%> construct a grid, which is usually used for a Shallow-Water type problem.\n
%> Edge information along with Dirichlet-type boundary conditions is read in 
%> from the fort.17-file. \n\n
%> Short overview of fort.17 files: \n
%> The first number is the total number of edges @f$numE@f$.\n
%> Then for the next @f$numE@f$ lines, i.e. for every edge the indices of the 
%> vertices at the end of the edge as well as the two elements that are 
%> separeted by this edge are given. For bounadry edges the second element index
%> must be zero.\n
%> Then for every element the three global edge indices are given in counter-
%> clockwise direction. \n
%> Next is the number of interior edges.\n
%> Then for every interior edge the global edge number is given.\n
%> Next is the number of land edges.\n
%> Then for every land edge the global edge number is given.\n
%> Next is the number of radiation edges.\n
%> Then for every radiation edge the global edge number is given.\n
%> Next is the number of river edges.\n
%> Then for every river edge the global edge number, the value of the free 
%> surface elevation and the normal and tangential velocities on this edge 
%> are given.\n
%> Next is the number of open sea edges.\n
%> Then for every open sea edge the global edge number is given.\n
%> For all frequencies and all open sea bounadry edges the amplitude and phase
%> in degrees are given.\n\n
%>
%> Since the grids in FESTUNG need to fulfill certain local indexing 
%> properties, such as that local edges of arbitrary index must be opposing 
%> local vertices of the same index, some of the information provided must be 
%> altered.
%>
%> @param  grid_file    The name of the fort.14-file.
%> @param  conn_file    The name of the fort.17-file.
%> @param  numForcingOS The number of forcing terms on the open sea boundary as
%>                      provided by the fort.15-parameter NBFR.
%> @param  isSpherical  <code>logical</code> scalar as provided by the fort.15-
%>                      parameter ICS that indicates if CPP projection is used 
%>                      to include curvature of the physical domain or if it is
%>                      approximately flat.
%> @param  projCenter   Projection center for CPP projection as provided by the 
%>                      fort.15-parameters SLAM0, SFEA0 in case curvature is 
%>                      considered. @f$[1 \times 2]@f$
%>
%> @retval g            The lists describing the geometric and topological 
%>                      properties of a triangulation.
%> @retval  depth       Positive values of the depth in each vertex of the grid.
%>                      @f$[numV \times 1]@f$
%> @retval  forcingOS   Amplitude and phase (in degrees) of each harmonic 
%>                      forcing function at each open sea boundary edge.
%>                      @f$[numForcingOS \times numEbdrOS \times 2]@f$
%> @retval  flowRateRiv Free surface elevation, normal and tangential velocity
%>                      components of river boundary conditions for each river 
%>                      edge. @f$[numEbdRiv \times 3]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Balthasar Reuter, Florian Frank, Vadym Aizinger
%>                      Modified 08/23/16 by Hennes Hajduk
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
function [g, depth, forcingOS, flowRateRiv, flowRate] = domainADCIRC(grid_file, conn_file, numForcingOS, isSpherical, projCenter)
% Set default values for input parameters
switch nargin
  case 2
    numForcingOS = 0;
    isSpherical = false;
    projCenter = [0 0];
  case 3
    isSpherical = false;
    projCenter = [0 0];
  case 4
    projCenter = [0 0];
end % switch
assert(exist(grid_file, 'file') == 2, ['Grid file "' grid_file '" does not exist!'])
assert(exist(conn_file, 'file') == 2, ['Connectivity file "' conn_file '" does not exist!'])
assert(isnumeric(numForcingOS) && round(numForcingOS) == numForcingOS && numForcingOS >= 0, 'Invalid number of open sea boundary forcings.')
validateattributes(projCenter, {'numeric'}, {'size', [1 2]}, mfilename, 'projectionCenter')

%% Grid file
% Open grid file and read name and data
fileID = fopen(grid_file, 'rt');
name = fgets(fileID);
data = cell2mat(textscan(fileID, '%f', 'CommentStyle', '!'));
fclose(fileID);

% Extract element and vertex counts
numT = data(1);
numV = data(2);
assert(round(numT) == numT && numT >= 1, 'Invalid number of elements.')
assert(round(numV) == numV && numV >= 3, 'Invalid number of vertices.')

% Extract coordinates
coordV = [ data(4:4:4*numV), data(5:4:4*numV+1) ];
if isSpherical
  % Transformation from spherical to cartesian coordinates
  r = 6378206.4;
  coordV(:,1) = r*pi/180 * (coordV(:,1) - projCenter(1)) * cos(pi/180 * projCenter(2));
  coordV(:,2) = r*pi/180 * coordV(:,2);
end % if
validateattributes(coordV, {'numeric'}, {'size', [numV 2]}, mfilename, 'coordV');

% Extract bathymetry
depth = data(6:4:4*numV+2);
validateattributes(depth, {'numeric'}, {'size', [numV 1], '>=', 0}, mfilename, 'depth');

% Extract connectivity
elements = reshape(data(4*numV+2+1 : 4*numV+2+5*numT), 5, numT).';
assert(isequal(elements(:,2), 3*ones(numT, 1)), 'Not supported, non-triangular elements found in input file.')
V0T = elements(:,3:5);
validateattributes(V0T, {'numeric'}, {'size', [numT 3], '>', 0, '<=', numV}, mfilename, 'V0T');

% Generate grid data
g = generateGridData(coordV, V0T);
g.name = name;

%% Connectivity file
% Open connectivity file and read data
fileID = fopen(conn_file, 'rt');
data = cell2mat(textscan(fileID, '%f', 'CommentStyle', '!'));

% Extract edge connectivity
numE = data(1);
assert(numE == g.numE, 'Invalid number of edges!')
V0E = [ data(3 : 5 : 5*numE-2), data(4 : 5 : 5*numE-1) ];
validateattributes(V0E, {'numeric'}, {'size', [numE 2], '>', 0, '<=', numV}, mfilename, 'V0E');
T0E = [ data(5 : 5 : 5*numE), data(6 : 5 : 5*numE+1) ];
validateattributes(T0E, {'numeric'}, {'size', [numE 2], '>=', 0, '<=', numT}, mfilename, 'T0E');
offset = 1 + 5 * numE + 1;
E0T = [ data(offset+1 : 4 : offset+4*numT-3), ...
        data(offset+2 : 4 : offset+4*numT-2), ...
        data(offset+3 : 4 : offset+4*numT-1) ];
validateattributes(E0T, {'numeric'}, {'size', [numT 3], '>', 0, '<=', numE}, mfilename, 'E0T');
offset = offset + 4 * numT;

% Interior edges
numEint = data(offset);
idxEint = data(offset+2 : 2 : offset+2*numEint);
validateattributes(idxEint, {'numeric'}, {'size', [numEint 1], '>', 0, '<=', numE}, mfilename, 'idxEint');
offset = offset + 1 + 2 * numEint;

% Land edges
numEbdrLand = data(offset);
idxEbdrLand = data(offset+2 : 2 : offset+2*numEbdrLand);
validateattributes(idxEbdrLand, {'numeric'}, {'size', [numEbdrLand 1], '>', 0, '<=', numE}, mfilename, 'idxEbdrLand');
offset = offset + 1 + 2 * numEbdrLand;

% Radiation edges
numEbdrRad = data(offset);
idxEbdrRad = data(offset+2 : 2 : offset+2*numEbdrRad);
validateattributes(idxEbdrRad, {'numeric'}, {'size', [numEbdrRad 1], '>', 0, '<=', numE}, mfilename, 'idxEbdrRad');
offset = offset + 1 + 2 * numEbdrRad;

% River edges
% This procedure is needed, since river edge boundary conditions are epecified in the form 
% J XI=ETRI(J). U=UNRI(J). V=UTRI(J)
numEbdrRiv = data(offset);
idxEbdrRiv = zeros(numEbdrRiv, 1);
flowRateRiv = zeros(numEbdrRiv, 3);
if numEbdrRiv > 0
  idxEbdrRiv(1) = data(offset+2);
  data = fgets(fileID);
  [xi, uv] = strtok(data);
  [u, v] = strtok(uv);
  flowRateRiv(1,:) = [ str2double(xi(4:end)), str2double(u(3:end)), str2double(v(4:end)) ];
  for j = 2 : numEbdrRiv
    data = cell2mat(textscan(fileID, '%f', 'CommentStyle', '!'));
    idxEbdrRiv(j) = data(2);
    data = fgets(fileID);
    [xi, uv] = strtok(data);
    [u, v] = strtok(uv);
    flowRateRiv(j,:) = [ str2double(xi(4:end)), str2double(u(3:end)), str2double(v(4:end)) ];
  end % for
  data = cell2mat(textscan(fileID, '%f', 'CommentStyle', '!'));
  offset = 1;
else
  offset = offset + 1;
end  %if
validateattributes(idxEbdrRiv, {'numeric'}, {'size', [numEbdrRiv 1], '>', 0, '<=', numE}, mfilename, 'idxEbdrRiv');

% Open Sea edges
numEbdrOS = data(offset);
idxEbdrOS = data(offset+2 : 2 : offset+2*numEbdrOS);
if numEbdrOS > 0
  validateattributes(idxEbdrOS, {'numeric'}, {'size', [numEbdrOS 1], '>', 0, '<=', numE}, mfilename, 'idxEbdrOS');
end % if
offset = offset + 1 + 2 * numEbdrOS;

% Open sea boundary forcings
forcingOS = zeros(numForcingOS, numEbdrOS, 2);
if numEbdrOS > 0
  for i = 1 : numForcingOS
    forcingOS(i, :, :) = [ data(offset : 2 : offset+2*numEbdrOS-2), data(offset+1 : 2 : offset+2*numEbdrOS-1) ];
    offset = offset + 2*numEbdrOS;
  end % for
end % if
% offset = offset + 2 * numEbdrOS;

% Flow edges
if length(data) >= offset
  numEbdrF = data(offset);
  idxEbdrF = data(offset+2 : 4 : offset+4*numEbdrF);
  flowRate = zeros(numEbdrF, 2);
  if numEbdrF > 0
    flowRate(:,1) = data(offset+3 : 4 : offset+4*numEbdrF);
    flowRate(:,2) = data(offset+4 : 4 : offset+4*numEbdrF);
    validateattributes(idxEbdrF, {'numeric'}, {'size', [numEbdrF 1], '>', 0, '<=', numE}, mfilename, 'idxEbdrF');
  end % if
else
  numEbdrF = 0;
  idxEbdrF = [];
  flowRate = zeros(0, 2);
end % if

% Close connectivity file
fclose(fileID);

%% Verify and store edge counts
assert( numEint + numEbdrLand + numEbdrRad + numEbdrRiv + numEbdrOS + numEbdrF == numE, 'Total number of edges does not match sum of edges of each type.')
markIntLand = intersect(idxEint, idxEbdrLand);
markIntRad = intersect(idxEint, idxEbdrRad);
markIntRiv = intersect(idxEint, idxEbdrRiv);
markIntOS = intersect(idxEint, idxEbdrOS);
markIntF = intersect(idxEint, idxEbdrF);
markLandRad = intersect(idxEbdrLand, idxEbdrRad);
markLandRiv = intersect(idxEbdrLand, idxEbdrRiv);
markLandOS = intersect(idxEbdrLand, idxEbdrOS);
markLandF = intersect(idxEbdrLand, idxEbdrF);
markRadRiv = intersect(idxEbdrRad, idxEbdrRiv);
markRadOS = intersect(idxEbdrRad, idxEbdrOS);
markRadF = intersect(idxEbdrRad, idxEbdrF);
markRivOS = intersect(idxEbdrRiv, idxEbdrOS);
markRivF = intersect(idxEbdrRiv, idxEbdrF);
markOSF = intersect(idxEbdrOS, idxEbdrF);
isUniquelyDefined = isempty(markIntLand) && isempty(markIntRad) && isempty(markIntRiv) && isempty(markIntOS) && ...
                    isempty(markIntF) && isempty(markLandRad) && isempty(markLandRiv) && isempty(markLandOS) && ...
                    isempty(markLandF) && isempty(markRadRiv) && isempty(markRadOS) && isempty(markRadF) && ...
                    isempty(markRivOS) && isempty(markRivF) && isempty(markOSF);
assert(isUniquelyDefined, 'At least one edge was given more than one type.')

% Store edge counts
g.numEint = numEint;
g.numEbdrL = numEbdrLand;
g.numEbdrRA = numEbdrRad;
g.numEbdrRI = numEbdrRiv;
g.numEbdrOS = numEbdrOS;
g.numEbdrF = numEbdrF;

%% Map edge indices and set edge types
% Sort edge-to-vertex mappings and search for identical rows
[markMember, adcirc2g] = ismember(sort(V0E,2), sort(g.V0E,2), 'rows');
validateattributes(markMember, {'logical'}, {'nonzero'}, mfilename, 'markMember')

% Set edge types
g.idE = NaN(g.numE, 1);
g.idE(adcirc2g(idxEint)) = 0;     % interior edges
g.idE(adcirc2g(idxEbdrLand)) = 1; % land boundary
g.idE(adcirc2g(idxEbdrRad)) = 2;  % radiation boundary
g.idE(adcirc2g(idxEbdrRiv)) = 3;  % river boundary
g.idE(adcirc2g(idxEbdrOS)) = 4;   % open sea boundary
g.idE(adcirc2g(idxEbdrF)) = 5;   % flow boundary
g.idE0T = g.idE(g.E0T);           % local edge IDs
assert(isequal(isnan(g.idE), zeros(numE, 1)), 'Missing edge type specification')

% Reorder boundary forcings
[~, I] = sort(adcirc2g(idxEbdrOS));
forcingOS = forcingOS(:,I,:);
[~, I] = sort(adcirc2g(idxEbdrRiv));
flowRateRiv = flowRateRiv(I,:);
[~, I] = sort(adcirc2g(idxEbdrF));
flowRate = flowRate(I,:);
end % function

