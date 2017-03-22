function problemData = solveStep(problemData, nStep)
% Obtain Runge-Kutta rule
[problemData.t, problemData.omega] = rungeKuttaSSP(problemData.ordRK, problemData.tau, (nStep - 1) * problemData.tau);

% Initialize solution vectors for RK steps
problemData.cDiscRK = cell(length(problemData.omega)+1, 3);
problemData.cDiscRK(1,:) = problemData.cDisc;

% Carry out RK steps
problemData.isSubSteppingFinished = false;
problemData = iterateSubSteps(problemData, nStep);

% Store solution
problemData.cDisc = problemData.cDiscRK(end, :);
end % function