function stepHandles = getStepHandles(problemName, stepList)
if nargin < 2
  [preprocessList, stepList, postprocessList, subStepList] = getStepLists();
  if isequal(cellfun(@(fun) exist([problemName filesep fun '.m'], 'file'), subStepList), 2 * ones(size(subStepList)))
    stepList = [ preprocessList(:); stepList(:); subStepList(:); postprocessList(:) ];
  else
    stepList = [ preprocessList(:); stepList(:); postprocessList(:) ];
  end % if
end % if
stepHandles = struct();
for nFunc = 1 : length(stepList)
  stepHandles.(stepList{nFunc}) = getFunctionHandle([problemName filesep stepList{nFunc}]);
end % for
% for nFunc = 1 : length(preprocessList)
%   stepHandles.(preprocessList{nFunc}) = getFunctionHandle([problemName filesep preprocessList{nFunc}]);
% end % for
% for nFunc = 1 : length(stepList)
%   stepHandles.(stepList{nFunc}) = getFunctionHandle([problemName filesep stepList{nFunc}]);
% end % for
% if isequal(cellfun(@(fun) exist([problemName filesep fun '.m'], 'file'), subStepList), 2 * ones(size(subStepList)))
%   for nFunc = 1 : length(subStepList)
%     stepHandles.(subStepList{nFunc}) = getFunctionHandle([problemName filesep subStepList{nFunc}]);
%   end % for
% end % if
% for nFunc = 1 : length(postprocessList)
%   stepHandles.(postprocessList{nFunc}) = getFunctionHandle([problemName filesep postprocessList{nFunc}]);
% end % for
end % function