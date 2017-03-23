function problemData = solveStep(problemData, nStep)
% Obtain Runge-Kutta rule
[problemData.t, problemData.omega] = rungeKuttaSSP(problemData.ordRK, problemData.tau, (nStep - 1) * problemData.tau);

% Carry out RK steps
problemData.isSubSteppingFinished = false;
problemData = iterateSubSteps(problemData, nStep);
end % function