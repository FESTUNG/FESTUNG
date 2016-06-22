function problemData = initializeProblem(problemData)
addpath('swe');
problemData.sweData = initializeProblem(problemData.sweData);
rmpath('swe');
addpath('transport');
problemData.transportData = initializeProblem(problemData.transportData);
rmpath('transport');
problemData.isFinished = problemData.sweData.isFinished || problemData.transportData.isFinished; % TODO ok? make consistent with postprocessStep
end % function