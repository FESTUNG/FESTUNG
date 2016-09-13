function visualizeDataSub( g , repLagr , varName, fileName, tLvl )

[K, N] = size(repLagr);
%% open file
fileName    = [fileName, '.', num2str(tLvl), '.vtu'];
file        = fopen(fileName, 'wt');
%% header
fprintf(file, '<?xml version="1.0"?>\n');
fprintf(file, '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">\n');
fprintf(file, '  <UnstructuredGrid>\n');
%% points & cells
switch N
    case {1, 4}
        P1 = reshape(g.coordV0Tsub(:,:,1)', 4*K, 1);
        P2 = reshape(g.coordV0Tsub(:,:,2)', 4*K, 1);
        numP = 4;
        id = 9;
    case 9
        P1 = reshape([g.coordV0Tsub(:,:,1), g.baryE0Tsub(:,[1,3,2,4],1)]', 8*K, 1);
        P2 = reshape([g.coordV0Tsub(:,:,2), g.baryE0Tsub(:,[1,3,2,4],2)]', 8*K, 1);
        numP = 8;
        id = 23;
        N = 8;
end  % switch
fprintf(file, '    <Piece NumberOfPoints="%d" NumberOfCells="%d">\n',K*numP,K);
fprintf(file, '      <Points>\n');
fprintf(file, '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">\n');
fprintf(file, '          %.3e %.3e %.3e\n',  [P1, P2, zeros(numP*K, 1)]');
fprintf(file, '        </DataArray>\n');
fprintf(file, '      </Points>\n');
fprintf(file, '      <Cells>\n');
fprintf(file, '        <DataArray type="Int32" Name="connectivity" format="ascii">\n');
fprintf(file,'           '); fprintf(file,'%d ', 0:K*numP-1);
fprintf(file, '\n        </DataArray>\n');
fprintf(file, '        <DataArray type="Int32" Name="offsets" format="ascii">\n');
fprintf(file,'           %d\n', numP:numP:numP*K);
fprintf(file, '        </DataArray>\n');
fprintf(file, '        <DataArray type="UInt8" Name="types" format="ascii">\n');
fprintf(file,'           %d\n', id*ones(K, 1));
fprintf(file, '        </DataArray>\n');
fprintf(file, '      </Cells>\n');
%% data
switch N
    case 1 % locally constant
        dataLagr = kron(repLagr, [1;1;1;1])';
    case 4 % locally quadratic
        dataLagr = reshape(repLagr', 1, K*N);
    case 8
        repLagr = repLagr(:,(1:N));
        dataLagr = reshape(repLagr', 1, K*N);
end  % switch
fprintf(file, '      <PointData Scalars="%s">\n', varName);
fprintf(file, '        <DataArray type="Float32" Name="%s" NumberOfComponents="1" format="ascii">\n', varName);
fprintf(file, '          %.3e\n', dataLagr);
fprintf(file, '        </DataArray>\n');
fprintf(file, '      </PointData>\n');
%% footer
fprintf(file, '    </Piece>\n');
fprintf(file, '  </UnstructuredGrid>\n');
fprintf(file, '</VTKFile>\n');
%% close file
fclose(file);
disp(['Data written to ' fileName])
end

