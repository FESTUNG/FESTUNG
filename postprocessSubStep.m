function problemData = postprocessSubStep(problemData, nStep, nSubStep)
% TODO add swe/postprocessingSubStep
addpath('transport');
problemData.transportData = postprocessSubStep(problemData.transportData, nStep, nSubStep);
rmpath('transport');
problemData.isSubSteppingFinished = problemData.transportData.isSubSteppingFinished; % TODO change if different RK schemes should be allowed
end % function
