function problemData = outputStep(problemData, nStep)
addpath('swe');
problemData.sweData = outputStep(problemData.sweData, nStep);
rmpath('swe');
addpath('transport');
problemData.transportData = outputStep(problemData.transportData, nStep);
rmpath('transport');
end % function

