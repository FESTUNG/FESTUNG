function problemData = solveSubStep(problemData, nStep, nSubStep)
problemData.sweData = problemData.swe_solveSubStep(problemData.sweData, nStep, nSubStep);
problemData.transportData = problemData.transport_solveSubStep(problemData.transportData, nStep, nSubStep);
end % function