function problemData = preprocessStep(problemData, nStep)
%% Apply mesh adaptation to free surface movement
problemData = problemData.fn_adaptMesh(problemData);
end % function