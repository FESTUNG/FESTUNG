function problemData = initializeProblem(problemData)
problemData.isFinished = false;
addpath('transport');
problemData.transportData = initializeProblem(problemData.transportData);
rmpath('transport');
end % function