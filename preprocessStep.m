function problemData = preprocessStep(problemData, nStep)
problemData.sweData = problemData.swe_preprocessStep(problemData.sweData, nStep);
problemData.transportData = problemData.transport_preprocessStep(problemData.transportData, nStep);
end % function
