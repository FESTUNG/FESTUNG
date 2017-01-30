function problemData = solveStep(problemData, nStep)
%% Convert representation matrix to representation vector
[barK,barN] = size(problemData.cDisc{1});
[K,N] = size(problemData.cDisc{2});
hSys = reshape(problemData.cDisc{1}.', barK*barN, 1);
cSys = cell(3,1);
cSys(2:3) = cellfun(@(c) reshape(c.', K*N, 1), problemData.cDisc(2:3), 'UniformOutput', false);
%% Solve for next time level
% Flux variables
qSys = cell(2,1);
for m = 1 : 2
  qSys{m} = problemData.globM \ ( -problemData.globJD{m} + ...
            (problemData.globH{m} - problemData.globQ{m}) * cSys{2} );
end % for m
% Water height
cSys{1} = hSys + problemData.tau * problemData.barGlobM \ ( problemData.globLh + problemData.barGlobJh - ...
            (problemData.barGlobP + problemData.barGlobPbdr) * cSys{2} + problemData.barGlobG * hSys );
% Horizontal velocity
cSys{2} = cSys{2} + problemData.tau * problemData.globM \ ( problemData.globLu + problemData.globJu + ...
            (problemData.globE{1} - problemData.globP{1} - problemData.globPbdr{1}) * cSys{2} + ...
            (problemData.globE{2} - problemData.globP{2} - problemData.globPbdr{2}) * cSys{3} + ...
            (problemData.globG{1} - problemData.globR{1} - problemData.globRbdr{1}) * qSys{1} + ...
            (problemData.globG{2} - problemData.globR{2} - problemData.globRbdr{2}) * qSys{2} + ...
            (problemData.gConst * problemData.tildeGlobH{1} - problemData.tildeGlobQ{1} - problemData.tildeGlobQbdr{1}) * hSys );
% Vertical velocity
cSys{3} = (problemData.globH{2} - problemData.globQup) \ (problemData.globJD{1} + problemData.globJh + ...
                          (-problemData.globH{1} + problemData.globQavg + problemData.tildeGlobP) * cSys{2} );
%% Convert representation vector to representation matrix
problemData.cDisc{1} = reshape(cSys{1}, barN, barK).';
problemData.cDisc(2:3) = cellfun(@(c) reshape(c, N, K).', cSys(2:3), 'UniformOutput', false);                        
end % function