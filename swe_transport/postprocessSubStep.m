function problemData = postprocessSubStep(problemData, nStep, nSubStep)
problemData.sweData = problemData.swe_postprocessSubStep(problemData.sweData, nStep, nSubStep);
problemData.transportData = problemData.transport_postprocessSubStep(problemData.transportData, nStep, nSubStep);

problemData.isSubSteppingFinished = problemData.transportData.isSubSteppingFinished; % TODO change if different RK schemes should be allowed
end % function
