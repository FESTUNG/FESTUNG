function problemData = solveStep(problemData, nStep)
KN = problemData.K * problemData.N;

% Obtain Runge-Kutta rule
[problemData.timeLvls, problemData.omega] = rungeKuttaSSP(problemData.ordRK, problemData.tau, (nStep - 1) * problemData.tau);

% Linearize solution vector
problemData.cDiscRK0 = cellfun(@(c) reshape(c.', [KN,1]), problemData.cDisc, 'UniformOutput', false);
problemData.cDiscRK = problemData.cDiscRK0;

% Initialize substepping
problemData.isSubSteppingFinished = false;
problemData = iterateSubSteps(problemData, nStep);
end % function