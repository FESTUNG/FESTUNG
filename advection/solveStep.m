function problemData = solveStep(problemData, nStep)
K = problemData.K;
N = problemData.N;

% Obtain Runge-Kutta rule
[problemData.t, problemData.omega] = rungeKuttaSSP(problemData.ordRK, problemData.tau, (nStep - 1) * problemData.tau);

% Initialize solution vectors for RK steps
problemData.cDiscRK = cell(length(problemData.omega)+1, 1); 
problemData.cDiscRK{1} = reshape(problemData.cDisc', [K*N 1]);

% Carry out RK steps
problemData.isSubSteppingFinished = false;
problemData = iterateSubSteps(problemData, nStep);
end % function