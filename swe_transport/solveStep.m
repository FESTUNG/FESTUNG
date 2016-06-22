function problemData = solveStep(problemData, nStep)
K = problemData.transportData.K;
N = problemData.transportData.N;

% Obtain Runge-Kutta rule % TODO: change to allow changes in time increment
[problemData.timeLvls, problemData.omega] = rungeKuttaSSP(problemData.sweData.schemeOrder, problemData.sweData.dt, ...
                                                          (nStep - 1) * problemData.sweData.dt); % TODO use this for SWE once substepping is implemented

problemData.transportData.timeLvls = problemData.timeLvls;
problemData.transportData.omega = problemData.omega;

problemData.transportData.cDiscRK = cell(length(problemData.omega)+1, problemData.transportData.numSpecies);
problemData.transportData.cDiscRK(1,:) = cellfun(@(c) reshape(c.', [K*N,1]), problemData.transportData.cDisc, 'UniformOutput', false);

problemData.isSubSteppingFinished = false;
problemData = iterateSubSteps(problemData, nStep);
end % function