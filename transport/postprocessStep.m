function problemData = postprocessStep(problemData, nStep)
K = problemData.K;
N = problemData.N;

%% Reshape and store solution in problemData
problemData.cDisc = cellfun(@(c) reshape(c, [N, K]).', problemData.cDiscRK(end,:), 'UniformOutput', false);
problemData.isFinished = nStep >= problemData.numSteps;
end % function