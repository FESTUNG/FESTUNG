function problemData = postprocessSubStep(problemData, nStep, nSubStep)
addpath('transport');
problemData.transportData = postprocessSubStep(problemData.transportData, nStep, nSubStep);
rmpath('transport');
problemData.isSubSteppingFinished = problemData.transportData.isSubSteppingFinished;
end % function
