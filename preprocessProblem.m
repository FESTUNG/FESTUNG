function problemData = preprocessProblem(problemData)
addpath('swe');
problemData.sweData = preprocessProblem(problemData.sweData);
rmpath('swe');

problemData.transportData.g = problemData.sweData.g;
problemData.transportData.K = problemData.sweData.K;
problemData.transportData.tau = problemData.sweData.dt;
problemData.transportData.velN = problemData.sweData.N;

addpath('transport');
problemData.transportData = preprocessProblem(problemData.transportData);
rmpath('transport');
end % function