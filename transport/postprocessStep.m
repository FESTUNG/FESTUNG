function problemData = postprocessStep(problemData, nStep)
%% Reshape and store solution in problemData
problemData.cDisc = cellfun(@(c) reshape(c, [problemData.N, problemData.K]).', problemData.cDiscRK, 'UniformOutput', false);
problemData.isFinished = nStep >= problemData.numSteps;
end % function