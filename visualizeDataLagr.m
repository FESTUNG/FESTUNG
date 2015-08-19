% Write output files in VTK- or Tecplot-format.
%
%===============================================================================
%> @file visualizeDataLagr.m
%>
%> @brief Write output files in VTK- or Tecplot-format.
%===============================================================================
%>
%> @brief Write output files in VTK- or Tecplot-format.
%>
%> Depending on the given <code>fileType</code> (<code>'vtk'</code> (default) or
%> <code>'tp'</code>), it writes a <code>.vtu</code> or <code>.plt</code> file
%> for the visualization of a discrete quantity in 
%> @f$\mathbb{P}_p(\mathcal{T}_h), p \in \{0, 1, 2\}@f$.
%>
%> The name of the generated file ist <code>fileName.tLvl.vtu</code> or
%> <code>fileName.tLvl.plt</code>, respectively, where <code>tLvl</code> stands
%> for time level.
%>
%> @note Although the VTK-format supports @f$p=2@f$, Paraview (4.2.0 at the time
%>       of writing) splits each triangle into four triangles and visualizes
%>       the function as piecewise linear.
%>
%> @note The Tecplot file format doesn't support higher order functions, 
%>       therefore each triangle is split into four triangles with linear
%>       representation.
%>
%> @par Example
%> @parblock
%> @code
%> g = generateGridData([0, -1; sqrt(3), 0; 0, 1; -sqrt(3), 0], [4,1,3; 1,2,3]);
%> g.idE = (abs(g.nuE(:,2)) > 0) .* ((g.nuE(:,1)>0) + (g.nuE(:,2)>0)*2+1);
%> fAlg = @(X1, X2) (X1<0).*(X1.^2 - X2.^2 - 1) + (X1>=0).*(-X1.^2 - X2.^2 + 1);
%> for N = [1, 3, 6]
%>   p = (sqrt(8*N+1)-3)/2;
%>   quadOrd = max(2*p, 1);
%>   computeBasesOnQuad(N);
%>   fDisc = projectFuncCont2DataDisc(g, fAlg, quadOrd, integrateRefElemPhiPhi(N));
%>   fLagr = projectDataDisc2DataLagr(fDisc);
%>   visualizeDataLagr(g, fLagr, 'funname', ['fDOF', int2str(N)], 1, 'vtk');
%> end
%> @endcode
%> produces the following output using Paraview:
%> @image html  visP0.png  "fDOF1.1.vtu with range [-2/3,1/3]" width=1cm
%> @image html  visP1.png  "fDOF3.1.vtu with range [-8/5,6/5]" width=1cm
%> @image html  visP2.png  "fDOF6.1.vtu with range [-2, 2]" width=1cm
%> @endparblock
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  dataLagr   The Lagrangian representation of the quantity, as produced
%>                    by <code>projectDataDisc2DataLagr</code> @f$[K \times N]@f$
%> @param  varName    The name of the quantity within the output file
%> @param  fileName   The basename of the output file
%> @param  tLvl       The time level
%> @param  fileType   (optional) The output format to be written 
%>                    (<code>'vtk'</code> (default) or <code>'tp'</code>).
%>
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
function visualizeDataLagr(g, dataLagr, varName, fileName, tLvl, fileType)
%% Deduce default arguments
if nargin < 6 || isempty(fileType)
  fileType = 'vtk';
end
%% Call correct function for writing file.
if strcmp(fileType, 'vtk')
  visualizeDataLagrVtk(g, dataLagr, varName, fileName, tLvl);
elseif strcmp(fileType, 'tp')
  visualizeDataLagrTp(g, dataLagr, varName, fileName, tLvl);
else
  error('Unknown file type: %s', fileType);
end % if
end % function
%
%> @brief Helper routine to write VTK-files.
function visualizeDataLagrVtk(g, dataLagr, varName, fileName, tLvl)
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
  case 3 % locally linear
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
%
%> @brief Helper routine to write Tecplot-files.
function visualizeDataLagrTp(g, dataLagr, varName, fileName, tLvl)
[K, N] = size(dataLagr);
%% Open file.
fileName = [fileName, '.', num2str(tLvl), '.plt'];
file     = fopen(fileName, 'wt'); % if this file exists, then overwrite
%% Header.
fprintf(file, 'TITLE="FESTUNG output file"\n');
fprintf(file, 'VARIABLES=X, Y, "%s"\n', varName);
%% Points and cells.
switch N
  case {1, 3}
    P1          = reshape(g.coordV0T(:, :, 1)', 3*K, 1);
    P2          = reshape(g.coordV0T(:, :, 2)', 3*K, 1);
    numT        = K;
    V0T         = 1:length(P1);
  case 6
    P1          = reshape([g.coordV0T(:,:,1), g.baryE0T(:,[3,1,2],1)]',6*K,1);
    P2          = reshape([g.coordV0T(:,:,2), g.baryE0T(:,[3,1,2],2)]',6*K,1);
    numT        = 4*K;
    V0T         = reshape([ 1:6:6*K; 4:6:6*K; 6:6:6*K;
                            2:6:6*K; 5:6:6*K; 4:6:6*K;
                            3:6:6*K; 6:6:6*K; 5:6:6*K;
                            4:6:6*K; 5:6:6*K; 6:6:6*K ], 3, numT);
end % switch
%% Data.
switch N
  case 1 % locally constant
    dataLagr = kron(dataLagr, [1;1;1])';
  case 3 % locally linear
    dataLagr = reshape(dataLagr', 1, K*N)';
  case 6 % locally quadratic (permutation of local edge indices due to TP format)
    dataLagr = reshape(dataLagr(:, [1,2,3,6,4,5])', 1, K*N)';
end % switch
%% Zone header.
fprintf(file, 'ZONE T="Time=%.3e", ', tLvl);
fprintf(file, 'N=%d, E=%d, ', length(P1), numT);
fprintf(file, 'ET=TRIANGLE, F=FEBLOCK, ');
fprintf(file, 'SOLUTIONTIME=%.3e\n\n', tLvl);
%% Point coordinates and data.
fprintf(file, '%.3e %.3e %.3e %.3e %.3e\n', P1);
fprintf(file, '\n\n');
fprintf(file, '%.3e %.3e %.3e %.3e %.3e\n', P2);
fprintf(file, '\n\n');
fprintf(file, '%.3e %.3e %.3e %.3e %.3e\n', dataLagr);
fprintf(file, '\n');
%% Connectivity.
fprintf(file, '%d %d %d\n', V0T);
%% Close file.
fclose(file);
disp(['Data written to ' fileName])
end % function
