function problemData = solveStep(problemData, nStep)

% TODO

cSys = cellfun(@(c) reshape(c.', [], 1), problemData.cDisc, 'UniformOutput', false);

cSys{3} = (problemData.globH{2} - problemData.globQup) \ (problemData.globJD{1} + problemData.globJH + ...
                          (-problemData.globH{1} + problemData.globQavg + problemData.tildeGlobP) * cSys{2} );

end % function