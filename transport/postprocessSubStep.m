function problemData = postprocessSubStep(problemData, ~, nSubStep)
problemData.isSubSteppingFinished = nSubStep >= length(problemData.omega);
end % function
