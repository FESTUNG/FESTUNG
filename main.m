function main(problemName)
%% Check given problem
validateattributes(problemName, {'char'},{'nonempty'}, mfilename, 'problemName')
assert(isdir(problemName), 'No directory for specified problem found.')
%% List of functions making up a problem description
preprocessList = { 'configureProblem'; 'preprocessProblem'; 'initializeProblem' };
stepList = { 'preprocessStep'; 'solveStep'; 'postprocessStep'; 'outputStep' };
postprocessList = { 'postprocessProblem' };
%% Check existence of all required functions
assert(isequal(cellfun(@(fun) exist([problemName '/' fun '.m'], 'file'), preprocessList), 2 * ones(size(preprocessList))), ...
  'Not all the required functions for the preprocessing of the problem found.')
assert(isequal(cellfun(@(fun) exist([problemName '/' fun '.m'], 'file'), stepList), 2 * ones(size(stepList))), ...
  'Not all the required functions for the problem steps found.')
assert(isequal(cellfun(@(fun) exist([problemName '/' fun '.m'], 'file'), postprocessList), 2 * ones(size(postprocessList))), ...
  'Not all the required functions for the postprocessing of the problem found.')
%% Start logging and time measurements
more off % Disable paging of output
tic % Start time measurement
diary([problemName '.log']) % Start logging
%% Add problem to search path
oldpath = addpath(problemName);
%% Execute problem
try
  %% Pre-process and initialize problem
  tPreprocess = tic;
  problemData = struct;
  for nFunc = 1 : length(preprocessList)
    problemData = feval(preprocessList{nFunc}, problemData);
  end % for
  fprintf('Pre-processing time: %g seconds.\n', toc(tPreprocess));
  %% Enter iterative loop
  tLoop = tic;
  for nStep = 1 : problemData.numSteps
    for nFunc = 1 : length(stepList)
      problemData = feval(stepList{nFunc}, problemData, nStep);
    end % for
  end % for
  fprintf('Loop time: %g seconds.\n', toc(tLoop));
  %% Post-process problem
  tPostprocess = tic;
  for nFunc = 1 : length(postprocessList)
    problemData = feval(postprocessList{nFunc}, problemData);
  end % for
  fprintf('Post-processing time: %g seconds.\n', toc(tPostprocess))
catch e
  %% End logging and restore original path
  diary off
  path(oldpath);
  rethrow(e)
end % try/catch
fprintf('Total computation time: %g seconds.\n', toc);
diary off
%% Restore original search path
path(oldpath);
end % function