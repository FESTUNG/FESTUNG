function problemData = solveStep(problemData, nStep)
%% Convert representation matrix to representation vector
[barK,barN] = size(problemData.cDisc{1});
[K,N] = size(problemData.cDisc{2});
cSys = cellfun(@(c) reshape(c.', [], 1), problemData.cDisc, 'UniformOutput', false);
%% Solve for next time level
% Flux variables
qSys = cell(2,1);
for m = 1 : 2
  qSys{m} = problemData.globM \ ( -problemData.globJu{m} + ...
            (problemData.globH{m} - problemData.globQ{m} - problemData.globQbdr{m}) * cSys{2} );
end % for m
% Water height
hSys = cSys{1} + problemData.tau * ( problemData.globLh + problemData.barGlobM \ (-problemData.barGlobJuh - problemData.barGlobKh + ...
            (problemData.barGlobG - problemData.barGlobP - problemData.barGlobPbdr) * cSys{1}) );
% Horizontal velocity
cSys{2} = cSys{2} + problemData.tau * ( problemData.globLu + problemData.globM \ (-problemData.globKu - ...
            problemData.globJh - problemData.globJuu - problemData.globJuw - problemData.globJq + ...
            (problemData.globE - problemData.globP - problemData.globPbdr) * cSys{2} + ...
            (problemData.globG{1} - problemData.globR{1} - problemData.globRbdr{1}) * qSys{1} + ...
            (problemData.globG{2} - problemData.globR{2} - problemData.globRbdr{2}) * qSys{2} + ...
            (problemData.tildeGlobH{1} - problemData.tildeGlobQ{1} - problemData.tildeGlobQbdr{1}) * cSys{1}) );
% Vertical velocity
cSys{3} = (problemData.globH{2} - problemData.globQup) \ (problemData.globJu{1} + problemData.globJw{2} + problemData.globKh + ...
                          (-problemData.globH{1} + problemData.globQavg + problemData.tildeGlobP + problemData.globQbdr{1}) * cSys{2} );
%% Convert representation vector to representation matrix
problemData.cDisc{1} = reshape(hSys, barN, barK).';
problemData.cDisc(2:3) = cellfun(@(c) reshape(c, N, K).', cSys(2:3), 'UniformOutput', false);                        
end % function