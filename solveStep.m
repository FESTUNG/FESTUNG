function problemData = solveStep(problemData, nStep)
% Obtain Runge-Kutta rule
[problemData.t, problemData.omega] = rungeKuttaExplicit(problemData.ordRK, problemData.tau, problemData.t0 + (nStep - 1) * problemData.tau);

% Carry out RK steps
problemData.isSubSteppingFinished = false;
problemData = iterateSubSteps(problemData, nStep, problemData.subStepHandles);
end % function