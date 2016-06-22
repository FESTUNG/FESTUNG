function problemData = preprocessStep(problemData, nStep)
addpath('swe');
problemData.sweData = preprocessStep(problemData.sweData, nStep);
rmpath('swe');
addpath('transport');
problemData.transportData = preprocessStep(problemData.transportData, nStep);
rmpath('transport');
end % function
