function problemData = outputStep(problemData, nStep)
problemData.sweData = problemData.swe_outputStep(problemData.sweData, nStep);
problemData.transportData = problemData.transport_outputStep(problemData.transportData, nStep);
end % function

