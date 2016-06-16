function problemData = preprocessStep(problemData, nStep)
% No preprocessing necessary.
addpath('transport');
problemData.transportData = preprocessStep(problemData.transportData, nStep);
rmpath('transport');
end % function
