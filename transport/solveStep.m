function problemData = solveStep(problemData, nStep)
K = problemData.K;
N = problemData.N;

% Obtain Runge-Kutta rule
[problemData.timeLvls, problemData.omega] = rungeKuttaSSP(problemData.ordRK, problemData.tau, (nStep - 1) * problemData.tau);

problemData.cDiscRK = cell(length(problemData.omega)+1, problemData.numSpecies);
problemData.cDiscRK(1,:) = cellfun(@(c) reshape(c.', [K*N,1]), problemData.cDisc, 'UniformOutput', false);

problemData.isSubSteppingFinished = false;
problemData = iterateSubSteps(problemData, nStep);
end % function