function problemData = postprocessStep(problemData, nStep)
problemData.sweData = problemData.swe_postprocessStep(problemData.sweData, nStep);
problemData.transportData = problemData.transport_postprocessStep(problemData.transportData, nStep);

problemData.isFinished = problemData.transportData.isFinished || problemData.transportData.isFinished;
end % function