% Generic main program for arbitrary problems. Call with problem name given
% as a function argument to start computation.

%===============================================================================
%> @file main.m
%>
%> @brief Generic main program for arbitrary problems. Call with problem 
%>        name given as a function argument to start computation.
%===============================================================================
%>
%> @brief Generic main program for arbitrary problems. Call with problem 
%>        name given as a function argument to start computation.
%>
%> This routine provides a generic interface for arbitrary problems.
%> It makes the assumption that any kind of problem can be subdivided into
%> a number of blocks that are executed after each other:
%> 
%>  1. Configuration: Specifies all relevant problem parameters (such as
%>     mesh width, number of time steps, boundary conditions, etc.).
%>     See template/configureProblem.m
%>  2. Pre-processing: Computes all required data fields that are not
%>     changing over the course of the computation. Typically, this
%>     contains the computational mesh, stationary matrix blocks, etc.
%>     See template/preprocessProblem.m
%>  3. Initialization: Fills data structures necessary for the main loop
%>     with the initial data necessary to start the computation.
%>     See template/initializeProblem.m
%>  4. Main-Loop: Iteratively compute a solution, e.g., at different
%>     time-levels. Can also consist of only a single iteration (for
%>     stationary problems). A new iteration is entered whenever the
%>     parameter <code>problemData.isFinished</code> is not set to
%>     <code>true</code>.
%>  5. Post-processing: Performs tasks after the computation of the final
%>     solution is done. Typically, this can be some kind of error 
%>     evaluation or similar tasks.
%>     See template/postprocessProblem.m
%>
%> The main loop itself is again subdivided into four parts:
%> 
%>  1. Pre-processing: Performs all tasks (e.g., assembly of matrix blocks,
%>     evaluation of boundary conditions, etc.) that are required for the
%>     computation of the next solution.
%>     See template/preprocessStep.m
%>  2. Solution: Computes the solution at the next step, e.g., at a new
%>     time level.
%>     See template/solveStep.m
%>  3. Post-processing: Performs any tasks that are necessary after
%>     computing the next solution, e.g., slope-limiting. This is usually
%>     also the place to decide, whether the main loop should be
%>     terminated. For this, <code>problemData.isFinished</code> must be
%>     set to <code>true</code>.
%>     See template/postprocessStep.m
%>  4. Output: Takes care of any tasks necessary for the visualization or
%>     interpretation of the computed solution, e.g., writing it to a file.
%>     See template/outputStep.m
%>
%> A struct <code>problemData</code> is passed to and returned from every
%> routine in the above steps, which allows to store problem data or 
%> computed values.
%>
%> To implement a new problem using this generic framework, simply create
%> a subfolder and implement all above mentioned routines in this folder.
%> Then call it as <code>main('folder-name')</code>.
%>
%> This function accepts either an optional struct parameter that contains 
%> configuration values to be used, or an optional sequence of keys and values
%> that are then used to initialize a struct.
%>
%> @par Example
%> @parblock
%> The following two uses are identical:
%> @code
%> problemData = struct();
%> problemData.p = 2;
%> problemData.tEnd = 1;
%> problemData = main('advection', problemData);
%> @endcode
%> @code
%> problemData = main('advection', 'p', 2, 'tEnd', 1);
%> @endcode
%> @endparblock
%>
%> @param  problemName  The name of the problem to be solved. A folder with
%>                      matching name must exist and provide all of the 
%>                      mentioned routines. @f$[\text{string}]@f$
%> @param  problemData  (optional) A struct containing configuration values
%>                      to be used by the solver. Thus, setdefault() should
%>                      be used in configureProblem() to honor given values.
%> @param  key          (optional) A string to identify the following value.
%> @param  value        (optional) The value corresponding to the preceding key.
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Balthasar Reuter, Florian Frank, Vadym Aizinger
%> 
%> @par License
%> @parblock
%> This program is free software: you can redistribute it and/or modify
%> it under the terms of the GNU General Public License as published by
%> the Free Software Foundation, either version 3 of the License, or
%> (at your option) any later version.
%>
%> This program is distributed in the hope that it will be useful,
%> but WITHOUT ANY WARRANTY; without even the implied warranty of
%> MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%> GNU General Public License for more details.
%>
%> You should have received a copy of the GNU General Public License
%> along with this program.  If not, see <http://www.gnu.org/licenses/>.
%> @endparblock
%
function varargout = main(problemName, varargin)
%% Check given problem
narginchk(1, inf)
nargoutchk(0, 1)
validateattributes(problemName, {'char'},{'nonempty'}, mfilename, 'problemName')
assert(isdir(problemName), 'No directory for specified problem found.')
if nargin == 2
  problemData = varargin{1};
  validateattributes(problemData, {'struct'}, {}, mfilename, 'problemData')
elseif nargin > 2
  assert(mod(nargin, 2) == 1, 'Optional arguments must take the form "key, value"')
  problemData = struct(varargin{:});
else
  problemData = struct;
end % if
problemData.problemName = problemName;
%% List of functions making up a problem description
[preprocessList, stepList, postprocessList] = getStepLists();
%% Check existence of all required functions
assert(isequal(cellfun(@(fun) exist([problemName filesep fun '.m'], 'file'), preprocessList), 2 * ones(size(preprocessList))), ...
  'Not all the required functions for the preprocessing of the problem found.')
assert(isequal(cellfun(@(fun) exist([problemName filesep fun '.m'], 'file'), stepList), 2 * ones(size(stepList))), ...
  'Not all the required functions for the problem steps found.')
assert(isequal(cellfun(@(fun) exist([problemName filesep fun '.m'], 'file'), postprocessList), 2 * ones(size(postprocessList))), ...
  'Not all the required functions for the postprocessing of the problem found.')
%% Start logging and time measurements, add problem to search path, and install exit handler
[tStartup, oldpath, cwd] = startupFestung(problemName);
cleanupObj = onCleanup(@() cleanupFestung(tStartup, oldpath, cwd));
%% Pre-process and initialize problem
tPreprocess = tic;
for nFunc = 1 : length(preprocessList)
  problemData = feval(preprocessList{nFunc}, problemData);
end % for
fprintf('Pre-processing time: %g seconds.\n', toc(tPreprocess));
%% Enter iterative loop
fprintf('Entering main loop.\n');
assert(isstruct(problemData) && isfield(problemData, 'isFinished') && islogical(problemData.isFinished), ...
  'Struct "problemData" must contain a logical variable "isFinished".');
tLoop = tic;
nStep = 0;
while ~problemData.isFinished
  nStep = nStep + 1;
  for nFunc = 1 : length(stepList)
    problemData = feval(stepList{nFunc}, problemData, nStep);
  end % for
end % while
tLoop = toc(tLoop);
fprintf('Loop time: %d iterations in %g seconds (on avg. %g seconds per iteration).\n', nStep, tLoop, tLoop/nStep);
%% Post-process problem
tPostprocess = tic;
for nFunc = 1 : length(postprocessList)
  problemData = feval(postprocessList{nFunc}, problemData);
end % for
fprintf('Post-processing time: %g seconds.\n', toc(tPostprocess))
%% Assign output variable
if nargout > 0
  varargout{1} = problemData;
end % if
end % function
%
function [tStartup, oldpath, cwd] = startupFestung(problemName)
more off % Disable paging of output
tStartup = tic; % Start time measurement
diaryName = [problemName '_' datestr(now, 'yyyymmdd-HHMMSS') '.log'];
diary(diaryName) % Start logging
fprintf( [ '\n' ...
'   __    __    __                      __    __    __ \n' ...
'  |  |  |  |  |  |                    |  |  |  |  |  |\n' ...
'  |  |__|  |__|  |                    |  |__|  |__|  |\n' ...
'  |              |   __    __    __   |              |\n' ...
'  |              |  |  |  |  |  |  |  |              |\n' ...
'  |              |__|  |__|  |__|  |__|              |\n' ...
'  |                                                  |\n' ...
'  |            Welcome to  F E S T U N G             |\n' ...
'  |                                                  |\n' ...
'  |      (c) 2014-%04s     Balthasar Reuter          |\n' ...
'  |                        Florian Frank             |\n' ...
'  |                        Vadym Aizinger            |\n' ...
'  |__________________________________________________|\n\n' ...
'   --------------------------------------------------\n' ...
'    Running problem "%s".\n' ...
'    Current date/time: %s.\n' ...
'    Logging output to "%s".\n' ...
'   --------------------------------------------------\n\n' ], ...
datestr(now,'yyyy'), problemName, datestr(now,'yyyy-mm-dd HH:MM:SS'), diaryName);
if isdeployed
  oldpath = path;
else
  oldpath = addpath([pwd filesep problemName], pwd);
end
cwd = pwd;
end % function
%
function cleanupFestung(tStartup, oldpath, cwd)
fprintf('Total computation time: %g seconds.\n', toc(tStartup));
diary off
path(oldpath);
cd(cwd);
end % function