function problemData = postprocessStep(problemData, nStep)
addpath('swe');
problemData.sweData = postprocessStep(problemData.sweData, nStep);
rmpath('swe');
addpath('transport');
problemData.transportData = postprocessStep(problemData.transportData, nStep);
rmpath('transport');
problemData.isFinished = problemData.transportData.isFinished || problemData.transportData.isFinished; % TODO ok, make consistent with initializeProblem
end % function

