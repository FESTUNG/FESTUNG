function [g, depth, forcingOS] = domainADCIRC(grid_file, conn_file, numForcingOS, isSpherical, projCenter)
% Set default values for input parameters
switch nargin
  case 3
    isSpherical = false;
  case 4
    projCenter = [0 0];
end % switch
assert(exist(grid_file, 'file') == 2, ['Grid file "' grid_file '" does not exist!'])
assert(exist(conn_file, 'file') == 2, ['Connectivity file "' conn_file '" does not exist!'])
assert(isscalar(numForcingOS) && round(numForcingOS) == numForcingOS && numForcingOS >= 0, 'Invalid number of open sea boundary forcings.')
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
E0T = [ data(offset+1 : 4 : offset+1+4*numT-3), ...
        data(offset+2 : 4 : offset+1+4*numT-2), ...
        data(offset+3 : 4 : offset+1+4*numT-1) ];
validateattributes(E0T, {'numeric'}, {'size', [numT 3], '>', 0, '<=', numE}, mfilename, 'E0T');
offset = offset + 4 * numT;

% Interior edges
numEint = data(offset);
idxEint = data(offset+2 : 2 : offset+1+2*numEint);
validateattributes(idxEint, {'numeric'}, {'size', [numEint 1], '>', 0, '<=', numE}, mfilename, 'idxEint');
offset = offset + 1 + 2 * numEint;

% Land edges
numEbdrLand = data(offset);
idxEbdrLand = data(offset+2 : 2 : offset+1+2*numEbdrLand);
validateattributes(idxEbdrLand, {'numeric'}, {'size', [numEbdrLand 1], '>', 0, '<=', numE}, mfilename, 'idxEbdrLand');
offset = offset + 1 + 2 * numEbdrLand;

% Radiation edges
numEbdrRad = data(offset);
idxEbdrRad = data(offset+2 : 2 : offset+1+2*numEbdrRad);
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
  flowRateRiv(1,:) = [ str2double(xi(4)), str2double(u(3)), str2double(v(4)) ];
  for j = 2 : numEbdrRiv
    data = cell2mat(textscan(fileID, '%f', 'CommentStyle', '!'));
    idxEbdrRiv(j) = data(2);
    data = fgets(fileID);
    [xi, uv] = strtok(data);
    [u, v] = strtok(uv);
    flowRateRiv(j,:) = [ str2double(xi(4)), str2double(u(3)), str2double(v(4)) ];
  end % for
  data = cell2mat(textscan(fileID, '%f', 'CommentStyle', '!'));
  offset = 1;
else
  offset = offset + 1;
end  %if
validateattributes(idxEbdrRiv, {'numeric'}, {'size', [numEbdrRiv 1], '>', 0, '<=', numE}, mfilename, 'idxEbdrRiv');

% Open Sea edges
numEbdrOS = data(offset);
idxEbdrOS = data(offset+2 : 2 : offset+1+2*numEbdrOS);
if numEbdrOS > 0
  validateattributes(idxEbdrOS, {'numeric'}, {'size', [numEbdrOS 1], '>', 0, '<=', numE}, mfilename, 'idxEbdrOS');
end % if
assert(numEbdrOS == 0 || numForcingOS > 0, 'Open sea boundary given but no boundary forcings.')
offset = offset + 1 + 2 * numEbdrOS;

% Open sea boundary forcings
forcingOS = zeros(numForcingOS, numEbdrOS, 2);
if numEbdrOS > 0
  for i = 1 : numForcingOS
    forcingOS(i, :, :) = [ data(offset : 2 : offset+2*numEbdrOS-2), data(offset+1 : 2 : offset+2*numEbdrOS-1) ];
    offset = offset + 2*numEbdrOS;
  end % for
end % if

% Close connectivity file
fclose(fileID);

%% Verify and store edge counts
assert( numEint + numEbdrLand + numEbdrRad + numEbdrRiv + numEbdrOS == numE, 'Total number of edges does not match sum of edges of each type.')
markIntLand = intersect(idxEint, idxEbdrLand);
markIntRad = intersect(idxEint, idxEbdrRad);
markIntRiv = intersect(idxEint, idxEbdrRiv);
markIntOS = intersect(idxEint, idxEbdrOS);
markLandRad = intersect(idxEbdrLand, idxEbdrRad);
markLandRiv = intersect(idxEbdrLand, idxEbdrRiv);
markLandOS = intersect(idxEbdrLand, idxEbdrOS);
markRadRiv = intersect(idxEbdrRad, idxEbdrRiv);
markRadOS = intersect(idxEbdrRad, idxEbdrOS);
markRivOS = intersect(idxEbdrRiv, idxEbdrOS);
isUniquelyDefined = isempty(markIntLand) && isempty(markIntRad) && isempty(markIntRiv) && isempty(markIntOS) && ...
                    isempty(markLandRad) && isempty(markLandRiv) && isempty(markLandOS) && ...
                    isempty(markRadRiv) && isempty(markRadOS) && isempty(markRivOS);
assert(isUniquelyDefined, 'At least one edge was given more than one type.')

% Store edge counts
g.numEint = numEint;
g.numEbdrLand = numEbdrLand;
g.numEbdrRad = numEbdrRad;
g.numEbdrRiv = numEbdrRiv;
g.numEbdrOS = numEbdrOS;

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
g.idE0T = g.idE(g.E0T);           % local edge IDs
assert(isequal(isnan(g.idE), zeros(numE, 1)), 'Missing edge type specification')
end % function

