function [preprocessList, stepList, postprocessList, subStepList] = getStepLists()
preprocessList = { 'configureProblem'; 'preprocessProblem'; 'initializeProblem' };
stepList = { 'preprocessStep'; 'solveStep'; 'postprocessStep'; 'outputStep' };
postprocessList = { 'postprocessProblem' };
subStepList = { 'preprocessSubStep'; 'solveSubStep'; 'postprocessSubStep' };
end % function