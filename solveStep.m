function problemData = solveStep(problemData, nStep)
K = problemData.transportData.K;
N = problemData.transportData.N;

% Obtain Runge-Kutta rule % TODO: change to allow changes in time increment
[problemData.timeLvls, problemData.omega] = rungeKuttaSSP(problemData.sweData.schemeOrder, problemData.sweData.dt, ...
                                                          (nStep - 1) * problemData.sweData.dt);

problemData.sweData.tLvls = problemData.timeLvls;
problemData.sweData.omega = problemData.omega;

problemData.sweData.cDiscRK0 = [ reshape(problemData.sweData.cDisc(:,:,1).', K*N, 1) ; ...
                                 reshape(problemData.sweData.cDisc(:,:,2).', K*N, 1) ; ...
                                 reshape(problemData.sweData.cDisc(:,:,3).', K*N, 1) ];
problemData.sweData.cDiscRK = problemData.sweData.cDiscRK0;

problemData.transportData.timeLvls = problemData.timeLvls;
problemData.transportData.omega = problemData.omega;

problemData.transportData.cDiscRK0 = cellfun(@(c) reshape(c.', [K*N,1]), problemData.transportData.cDisc, 'UniformOutput', false);
problemData.transportData.cDiscRK = problemData.transportData.cDiscRK0;

problemData.isSubSteppingFinished = false;
problemData = iterateSubSteps(problemData, nStep);
end % function