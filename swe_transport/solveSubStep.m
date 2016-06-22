function problemData = solveSubStep(problemData, nStep, nSubStep)
% TODO swe substepping once implemented
addpath('transport');
problemData.transportData = solveSubStep(problemData.transportData, nStep, nSubStep);
rmpath('transport');
end % function
