function problemData = outputStep(problemData, nStep)
% No output necessary.
addpath('transport');
problemData.transportData = outputStep(problemData.transportData, nStep);
rmpath('transport');
end % function

