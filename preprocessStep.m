function problemData = preprocessStep(problemData, nStep)
% Apply mesh adaptation to free surface movement
problemData = problemData.fn_adaptFreeSurface(problemData);

% Copy last time step
problemData.cDiscRK(1, :) = problemData.cDiscRK(end, :);
end % function