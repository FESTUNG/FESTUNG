function problemData = postprocessProblem(problemData)
addpath('swe');
problemData.sweData = postprocessProblem(problemData.sweData);
rmpath('swe');
addpath('transport');
problemData.transportData = postprocessProblem(problemData.transportData);
rmpath('transport');
end % function

