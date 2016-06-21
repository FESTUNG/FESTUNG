function problemData = solveStep(problemData, nStep)
K = problemData.K;
dataN = problemData.transportData.N;

% Obtain Runge-Kutta rule
[problemData.timeLvls, problemData.omega] = rungeKuttaSSP(problemData.ordRK, problemData.tau, (nStep - 1) * problemData.tau);

problemData.transportData.timeLvls = problemData.timeLvls;
problemData.transportData.omega = problemData.omega;

problemData.transportData.cDiscRK = cell(length(problemData.omega)+1, problemData.transportData.numSpecies);
problemData.transportData.cDiscRK(1,:) = cellfun(@(c) reshape(c.', [K*dataN,1]), problemData.transportData.cDisc, 'UniformOutput', false);

problemData.isSubSteppingFinished = false;
problemData = iterateSubSteps(problemData, nStep);
end % function