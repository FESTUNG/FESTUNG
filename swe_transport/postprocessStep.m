function problemData = postprocessStep(problemData, nStep)
addpath('transport');
problemData.transportData = postprocessStep(problemData.transportData, nStep);
rmpath('transport');
problemData.isFinished = problemData.transportData.isFinished;
end % function

