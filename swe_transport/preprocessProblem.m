function problemData = preprocessProblem(problemData)
problemData.sweData = execin('swe/preprocessProblem', problemData.sweData);

problemData.transportData.g = problemData.sweData.g;
problemData.transportData.K = problemData.sweData.K;
problemData.transportData.tau = problemData.sweData.dt;
problemData.transportData.velN = problemData.sweData.N;

problemData.transportData = execin('transport/preprocessProblem', problemData.transportData);

% only created function handles for routines that are called repeatedly
problemData.swe_preprocessStep = getFunctionHandle('swe/preprocessStep');
problemData.swe_solveStep = getFunctionHandle('swe/solveStep');
problemData.swe_preprocessSubStep = getFunctionHandle('swe/preprocessSubStep');
problemData.swe_solveSubStep = getFunctionHandle('swe/solveSubStep');
problemData.swe_postprocessSubStep = getFunctionHandle('swe/postprocessSubStep');
problemData.swe_postprocessStep = getFunctionHandle('swe/postprocessStep');
problemData.swe_outputStep = getFunctionHandle('swe/outputStep');
problemData.transport_preprocessStep = getFunctionHandle('transport/preprocessStep');
problemData.transport_solveStep = getFunctionHandle('transport/solveStep');
problemData.transport_preprocessSubStep = getFunctionHandle('transport/preprocessSubStep');
problemData.transport_solveSubStep = getFunctionHandle('transport/solveSubStep');
problemData.transport_postprocessSubStep = getFunctionHandle('transport/postprocessSubStep');
problemData.transport_postprocessStep = getFunctionHandle('transport/postprocessStep');
problemData.transport_outputStep = getFunctionHandle('transport/outputStep');
end % function