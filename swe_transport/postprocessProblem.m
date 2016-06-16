function problemData = postprocessProblem(problemData)
% No postprocessing necessary.
addpath('transport');
problemData.transportData = postprocessProblem(problemData.transportData);
rmpath('transport');
end % function

