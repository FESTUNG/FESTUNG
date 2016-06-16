function problemData = solveSubStep(problemData, nStep, nSubStep)

addpath('transport');
problemData.transportData = solveSubStep(problemData.transportData, nStep, nSubStep);
rmpath('transport');
end % function
