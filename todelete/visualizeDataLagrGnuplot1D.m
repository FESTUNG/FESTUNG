%
function visualizeDataLagrGnuplot1D(g, dataLagr, varName, fileName, tLvl, fileTypes, vecName)
    visualizeDataLagrGnuplot(g, dataLagr);
end % function
%

function visualizeDataLagrGnuplot(g, dataLagr)
%% Open file.
fileName = ['test', '.plt'];
file     = fopen(fileName, 'wt'); % if this file exists, then overwrite

len = size(dataLagr, 2);
Q1 = [1:-1/(len-1):0];

%% TODO Make this loop disappear
for i = 1:g.numE
    localIdx = g.E0E(i,1);
    adjTri = g.T0E(i, 1);
    
    F1 = @(X1, X2) g.B(adjTri,1,1)*X1 + g.B(adjTri,1,2)*X2 + g.coordV0T(adjTri,1,1)*ones(size(X1));
    F2 = @(X1, X2) g.B(adjTri,2,1)*X1 + g.B(adjTri,2,2)*X2 + g.coordV0T(adjTri,1,2)*ones(size(X1));

    [x1, x2] = gammaMap( localIdx, Q1 );
    
    X1 = F1(x1, x2);
    X2 = F2(x1, x2);
    
    for j=1:max(size(Q1))
        fprintf(file, '%.12e %.12e %.12e\n',  [ X1(j), X2(j), dataLagr(i, j)]');
    end    
end % for

%% Close file.
fclose(file);
disp(['Data written to ' fileName])
end
