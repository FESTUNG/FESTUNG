function problemData = postprocessStep(problemData, nStep)
K = problemData.K;
N = problemData.N;

%% Reshape and store solution in problemData
problemData.cDisc = reshape(problemData.cDiscRK{end}, N, K)';
problemData.isFinished = nStep >= problemData.numSteps;
end % function

