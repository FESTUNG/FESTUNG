% Write output files in VTK- or Tecplot-format.

%===============================================================================
%> @file
%>
%> @brief Write output files in VTK- or Tecplot-format.
%===============================================================================
%>
%> @brief Write output files in VTK- or Tecplot-format.
%>
%> Depending on the given <code>fileTypes</code> (<code>'vtk'</code> (default) or
%> <code>'tec'</code>), it writes a <code>.vtu</code> or <code>.plt</code> file
%> for the visualization of a discrete quantity in 
%> @f$\mathbb{P}_p(\mathcal{T}_h), p \in \{0, 1, 2\}@f$.
%> Multiple file types can be given at the same time, see the example
%> below.
%>
%> Multiple data sets can be written to a single file. To do so, specify
%> the different data sets as a cell array. Additionally, multiple data
%> sets can be grouped into a vector and written as such (VTK only).
%>
%> The name of the generated file ist <code>fileName.tLvl.vtu</code> or
%> <code>fileName.tLvl.plt</code>, respectively, where <code>tLvl</code> stands
%> for time level.

%> @par Example
%> @parblock
%> Assume a grid, <code>g</code>, and a cell array with discontinuous data,
%> <code>cDisc</code>, exist already.
%>
%> Then, the following code writes a single output file with a scalar field
%> for 'h' and a vector field 'velocity' with components 'u' and 'v'.
%> @code
%> varName = { 'h', 'u', 'v' };
%> vecName = struct('velocity', {{'u','v'}});
%> dataLagr = { projectDataDisc2DataLagrTensorProduct(cDisc{1}), ...
%>              projectDataDisc2DataLagrTensorProduct(cDisc{2}), ...
%>              projectDataDisc2DataLagrTensorProduct(cDisc{3}) };
%> visualizeDataLagrTetra(g, dataLagr, varName, 'solution', 1, {'vtk', 'tec'}), vecName);
%> end
%> @endcode
%> @endparblock
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  dataLagr   The Lagrangian representation of the quantity, as produced
%>                    by <code>projectDataDisc2DataLagr</code> @f$[K \times N]@f$.
%>                    Can be a cell array with multiple data sets.
%> @param  varName    The name of the quantity within the output file
%>                    If multiple data sets are given, this must be a cell
%>                    array with the same dimension as dataLagr, holding
%>                    the names of each data set.
%> @param  fileName   The basename of the output file
%> @param  tLvl       The time level
%> @param  fileTypes  (optional) The output format to be written 
%>                    (<code>'vtk'</code> (default) or <code>'tec'</code>).
%> @param  vecName    (optional) A struct that allows to define vectorial
%>                    output (for VTK only). Field names provide the vector
%>                    names and field values are cell arrays with variable
%>                    names.
%>
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
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
function visualizeDataLagrTetra(g, dataLagr, varName, fileName, tLvl, fileTypes, vecName)
%% Deduce default arguments
if nargin < 6 || isempty(fileTypes)
  fileTypes = 'vtk';
end % if
if nargin < 7
  vecName = struct;
end % if
if ~iscell(dataLagr)
  dataLagr = { dataLagr };
end % if
if ~iscell(varName)
  varName = { varName };
end % if
if ~iscell(fileTypes)
  fileTypes = cellstr(fileTypes);
end % if
if size(fileTypes,1) > size(fileTypes,2)
  fileTypes = transpose(fileTypes);
end % if
%% Check function arguments
assert(isequal(size(dataLagr), size(varName)), 'Number of data sets in dataLagr does not match varName')
assert(all(cellfun(@(c) size(c, 1), dataLagr) == g.numT), 'Wrong number of elements in dataLagr')
assert(all(cellfun(@(c) size(c, 2), dataLagr) == size(dataLagr{1}, 2)), 'All data sets in dataLagr must have same approximation order')
%% Ensure target directory exists
[dirName,~,~] = fileparts(fileName);
if ~isempty(dirName) && ~isdir(dirName)
  mkdir(dirName);
end % if
%% Call correct function for writing file.
for fileType = fileTypes
  if strcmp(fileType, 'vtk')
    visualizeDataLagrVtk(g, dataLagr, varName, fileName, tLvl, vecName);
  elseif strcmp(fileType, 'tec')
    visualizeDataLagrTec(g, dataLagr, varName, fileName, tLvl);
  else
    error('Unknown file type: %s', fileType);
  end % if
end
end % function
%
%> @brief Helper routine to write VTK-files.
function visualizeDataLagrVtk(g, dataLagr, varName, fileName, tLvl, vecName)
[K, N] = size(dataLagr{1});
vecNames = fieldnames(vecName);
isVec = ~isempty(vecNames);
maskVec = false(size(varName));
for i = 1 : length(vecNames)
  maskVec = maskVec | ismember(varName, vecName.(vecNames{i}));
end % for
isScalar = ~all(maskVec);
%% Open file.
fileName = [fileName, '.', num2str(tLvl), '.vtu'];
file     = fopen(fileName, 'wt'); % if this file exists, then overwrite
%% Header.
fprintf(file, '<?xml version="1.0"?>\n');
fprintf(file, '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">\n');
fprintf(file, '  <UnstructuredGrid>\n');
%% Points and cells.
switch N
  case {1, 4}
    P1          = reshape(g.coordV0T(:, :, 1)', 4*K, 1);
    P2          = reshape(g.coordV0T(:, :, 2)', 4*K, 1);
    numP        = 4; % number of local points
    id          = 9; % vtk ID for linear polynomials
  case 9
    P1          = reshape([g.coordV0T(:,:,1), g.baryE0T(:,[1,3,2,4],1)]',8*K,1);
    P2          = reshape([g.coordV0T(:,:,2), g.baryE0T(:,[1,3,2,4],2)]',8*K,1);
    numP        = 8; % number of local points
    id          = 23; % vtk ID for quadratic polynomials
end % switch
fprintf(file, '    <Piece NumberOfPoints="%d" NumberOfCells="%d">\n',K*numP,K);
fprintf(file, '      <Points>\n');
fprintf(file, '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">\n');
fprintf(file, '          %.12e %.12e %.1e\n',  [P1, P2, zeros(numP*K, 1)]');
fprintf(file, '        </DataArray>\n');
fprintf(file, '      </Points>\n');
fprintf(file, '      <Cells>\n');
fprintf(file, '        <DataArray type="Int32" Name="connectivity" format="ascii">\n');
fprintf(file, '           '); fprintf(file,'%d ', 0:K*numP-1);
fprintf(file, '\n        </DataArray>\n');
fprintf(file, '        <DataArray type="Int32" Name="offsets" format="ascii">\n');
fprintf(file, '           %d\n', numP:numP:numP*K);
fprintf(file, '        </DataArray>\n');
fprintf(file, '        <DataArray type="UInt8" Name="types" format="ascii">\n');
fprintf(file, '           %d\n', id*ones(K, 1));
fprintf(file, '        </DataArray>\n');
fprintf(file, '      </Cells>\n');
%% Data.
if isVec && isScalar
  fprintf(file, '      <PointData Scalars="%s" Vectors="%s">\n', varName{1}, vecNames{1});
elseif isScalar
  fprintf(file, '      <PointData Scalars="%s">\n', varName{1});
else
  fprintf(file, '      <PointData Vectors="%s">\n', vecNames{1});
end %if
for i = 1 : numel(dataLagr)
  switch N
    case 1 % locally constant
      dataLagr{i} = kron(dataLagr{i}, [1;1;1;1])';
    case 4 % locally linear
      dataLagr{i} = reshape(dataLagr{i}', 1, K*N);
    case 9 % locally quadratic
      dataLagr{i} = reshape(dataLagr{i}(:, 1:numP)', 1, K*numP);
  end % switch
end % for
for i = 1 : length(vecNames)
  maskComp = ismember(varName, vecName.(vecNames{i}));
  numComp = sum(maskComp);
  fprintf(file, '        <DataArray type="Float32" Name="%s" NumberOfComponents="3" format="ascii">\n', vecNames{i});
  fprintf(file, '          %.9e %.9e %.9e\n', [cell2mat(reshape(dataLagr(maskComp), [], 1)).', zeros(length(dataLagr{1}), 3-numComp)]');
  fprintf(file, '        </DataArray>\n');
end % for
for i = find(~maskVec)
  fprintf(file, '        <DataArray type="Float32" Name="%s" NumberOfComponents="1" format="ascii">\n', varName{i});
  fprintf(file, '          %.9e\n', dataLagr{i});
  fprintf(file, '        </DataArray>\n');
end % for
fprintf(file, '      </PointData>\n');
%% Footer.
fprintf(file, '    </Piece>\n');
fprintf(file, '  </UnstructuredGrid>\n');
fprintf(file, '</VTKFile>\n');
%% Close file.
fclose(file);
disp(['Data written to ' fileName])
end % function
%
%> @brief Helper routine to write Tecplot-files.
function visualizeDataLagrTec(g, dataLagr, varName, fileName, tLvl)
[K, N] = size(dataLagr{1});
%% Open file.
fileName = [fileName, '.', num2str(tLvl), '.plt'];
file     = fopen(fileName, 'wt'); % if this file exists, then overwrite
%% Header.
fprintf(file, 'TITLE="FESTUNG output file"\n');
fprintf(file, ['VARIABLES=X, Y', repmat(', "%s"', 1, numel(varName)), '\n'], varName{:});
%% Points and cells.
switch N
  case {1, 4}
    P1          = reshape(g.coordV0T(:, :, 1)', 4*K, 1);
    P2          = reshape(g.coordV0T(:, :, 2)', 4*K, 1);
    numT        = K;
    V0T         = 1:length(P1);
  case 9
    P1          = reshape([g.coordV0T(:,:,1), g.baryE0T(:,[1,3,2,4],1), sum(g.coordV0T(:,:,1), 2)/4]',9*K,1);
    P2          = reshape([g.coordV0T(:,:,2), g.baryE0T(:,[1,3,2,4],2), sum(g.coordV0T(:,:,2), 2)/4]',9*K,1);
    numT        = 4*K;
    V0T         = reshape([ 1:9:9*K; 5:9:9*K; 9:9:9*K; 8:9:9*K;
                            5:9:9*K; 2:9:9*K; 6:9:9*K; 9:9:9*K;
                            9:9:9*K; 6:9:9*K; 3:9:9*K; 7:9:9*K;
                            8:9:9*K; 9:9:9*K; 7:9:9*K; 4:9:9*K ], 4, numT);
end % switch
%% Data.
for i = 1 : numel(dataLagr)
  switch N
    case 1 % locally constant
      dataLagr{i} = kron(dataLagr{i}, [1;1;1;1])';
    case 4 % locally linear
      dataLagr{i} = reshape(dataLagr{i}', 1, K*N);
    case 9 % locally quadratic
      dataLagr{i} = reshape(dataLagr{i}', 1, K*N);
  end % switch
end % for
%% Zone header.
fprintf(file, 'ZONE T="Time=%.3e", ', tLvl);
fprintf(file, 'N=%d, E=%d, ', length(P1), numT);
fprintf(file, 'ET=QUADRILATERAL, F=FEBLOCK, ');
fprintf(file, 'SOLUTIONTIME=%.3e\n\n', tLvl);
%% Point coordinates and data.
fprintf(file, '%.12e %.12e %.12e %.12e %.12e\n', P1);
fprintf(file, '\n\n');
fprintf(file, '%.12e %.12e %.12e %.12e %.12e\n', P2);
fprintf(file, '\n\n');
fprintf(file, '%.9e %.9e %.9e %.9e %.9e\n', dataLagr{:});
fprintf(file, '\n');
%% Connectivity.
fprintf(file, '%d %d %d\n', V0T);
%% Close file.
fclose(file);
disp(['Data written to ' fileName])
end % function