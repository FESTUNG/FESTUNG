% This file is part of FESTUNG 
% Copyright (C) 2014 Florian Frank, Balthasar Reuter, Vadym Aizinger
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
function visualizeDataLagr(g, dataLagr, varName, fileName, tLvl)
[K, N] = size(dataLagr);
%% Open file.
fileName = [fileName, '.', num2str(tLvl), '.vtu'];
file     = fopen(fileName, 'wt'); % if this file exists, then overwrite
%% Header.
fprintf(file, '<?xml version="1.0"?>\n');
fprintf(file, '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">\n');
fprintf(file, '  <UnstructuredGrid>\n');
%% Points and cells.
switch N
  case {1, 3}
    P1          = reshape(g.coordV0T(:, :, 1)', 3*K, 1);
    P2          = reshape(g.coordV0T(:, :, 2)', 3*K, 1);
    numP        = 3; % number of local points
    id          = 5; % vtk ID for linear polynomials
  case 6
    P1          = reshape([g.coordV0T(:,:,1), g.baryE0T(:,[3,1,2],1)]',6*K,1);
    P2          = reshape([g.coordV0T(:,:,2), g.baryE0T(:,[3,1,2],2)]',6*K,1);
    numP        = 6; % number of local points
    id          = 22; % vtk ID for quadratic polynomials
end % switch
fprintf(file, '    <Piece NumberOfPoints="%d" NumberOfCells="%d">\n',K*numP,K);
fprintf(file, '      <Points>\n');
fprintf(file, '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">\n');
fprintf(file, '          %.3e %.3e %.3e\n',  [P1, P2, zeros(numP*K, 1)]');
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
switch N
  case 1 % locally constant
    dataLagr = kron(dataLagr, [1;1;1])';
  case 3 % locally quadratic
    dataLagr = reshape(dataLagr', 1, K*N);
  case 6 % locally quadratic (permutation of local edge indices due to vtk format)
    dataLagr = reshape(dataLagr(:, [1,2,3,6,4,5])', 1, K*N);
end % switch
fprintf(file, '      <PointData Scalars="%s">\n', varName);
fprintf(file, '        <DataArray type="Float32" Name="%s" NumberOfComponents="1" format="ascii">\n', varName);
fprintf(file, '          %.3e\n', dataLagr);
fprintf(file, '        </DataArray>\n');
fprintf(file, '      </PointData>\n');
%% Footer.
fprintf(file, '    </Piece>\n');
fprintf(file, '  </UnstructuredGrid>\n');
fprintf(file, '</VTKFile>\n');
%% Close file.
fclose(file);
disp(['Data written to ' fileName])
end % function
