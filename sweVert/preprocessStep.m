function problemData = preprocessStep(problemData, nStep)
problemData.cDiscRK(1, :) = problemData.cDiscRK(end, :);

% Apply mesh adaptation to free surface movement
problemData = problemData.fn_adaptMesh(problemData);
end % function