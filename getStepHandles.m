function stepHandles = getStepHandles(problemName)
[preprocessList, stepList, postprocessList, subStepList] = getStepLists();
stepHandles = struct();
for nFunc = 1 : length(preprocessList)
  stepHandles.(preprocessList{nFunc}) = getFunctionHandle([problemName filesep preprocessList{nFunc}]);
end % for
for nFunc = 1 : length(stepList)
  stepHandles.(stepList{nFunc}) = getFunctionHandle([problemName filesep stepList{nFunc}]);
end % for
if isequal(cellfun(@(fun) exist([problemName filesep fun '.m'], 'file'), subStepList), 2 * ones(size(subStepList)))
  for nFunc = 1 : length(subStepList)
    stepHandles.(subStepList{nFunc}) = getFunctionHandle([problemName filesep subStepList{nFunc}]);
  end % for
end % if
for nFunc = 1 : length(postprocessList)
  stepHandles.(postprocessList{nFunc}) = getFunctionHandle([problemName filesep postprocessList{nFunc}]);
end % for
end % function