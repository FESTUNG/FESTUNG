function problemData = preprocessStep(problemData, nStep)
%% Apply mesh adaptation to free surface movement
problemData = adaptMesh(problemData);
end % function