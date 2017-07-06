% Fills the problemData-struct with all basic configuration options.
% Problem parameters are to be modified inside this routine.

%===============================================================================
%> @file darcyVert_sweVert/configureProblem.m
%>
%> @brief Fills the problemData-struct with all basic configuration options.
%>        Problem parameters are to be modified inside this routine.
%===============================================================================
%>
%> @brief Fills the problemData-struct with all basic configuration options.
%>        Problem parameters are to be modified inside this routine.
%>
%> This routine is called before any other function for the problem.
%> It should define all problem parameters.
%>
%> The only requirement for this struct is to provide a boolean parameter 
%> <code>problemData.isFinished</code> that specifies whether the 
%> iterative solver has finished (i.e., the parameters value is 
%> <code>true</code>).
%> Typically, this parameter is introduced in initializeProblem() and
%> updated in postprocessStep() according to the progress of the solver.
%>
%> See main() for more details about the solver structure.
%>
%> @param  problemData  A (probably) empty struct with problem parameters.
%>                      @f$[\text{struct}]@f$
%>
%> @retval problemData  A struct with all necessary parameters and definitions
%>                      for the problem description. @f$[\text{struct}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
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
function problemData = configureProblem(problemData)
%% Parameters.
% Name of testcase
problemData = setdefault(problemData, 'testcase', 'coupling');

% Enable coupling
problemData = setdefault(problemData, 'isCouplingDarcy', true);
problemData = setdefault(problemData, 'isCouplingSWE', true);

% Number of elements in x- and y-direction
problemData = setdefault(problemData, 'numElem', [128, 64]);

% Local polynomial approximation order (0 to 5)
problemData = setdefault(problemData, 'p', 1);

% Order of quadrature rule
problemData = setdefault(problemData, 'qOrd', 2*problemData.p + 1);

% Time stepping parameters
problemData = setdefault(problemData, 't0', 0);  % start time
problemData = setdefault(problemData, 'tEnd', 0.01);  % end time
problemData = setdefault(problemData, 'numSteps', 50);  % number of time steps
problemData = setdefault(problemData, 'numSubSteps', 10); % number of free-flow steps per sub-surface step

% Visualization settings
problemData = setdefault(problemData, 'isVisGrid', false);  % visualization of grid
problemData = setdefault(problemData, 'isVisSol', false);  % visualization of solution
problemData = setdefault(problemData, 'outputFrequency', 100); % no visualization of every timestep

couplingString = '_';
if problemData.isCouplingDarcy, couplingString = [couplingString 'couplingDarcy']; end
if problemData.isCouplingSWE, couplingString = [couplingString 'couplingSWE']; end

problemData = setdefault(problemData, 'outputBasename', [problemData.testcase couplingString]); 
problemData = setdefault(problemData, 'outputTypes', { 'vtk' });  % Type of visualization files ('vtk, 'tec')

%% Function handles for steps of the sub-problems
problemData.darcySteps = getStepHandles('darcyVert');
problemData.sweSteps = getStepHandles('sweVert');

%% Sub-surface problem
problemData.darcyData = struct;
problemData.darcyData.problemName = 'darcyVert';
problemData.darcyData.testcase = problemData.testcase;

problemData.darcyData.numElem = problemData.numElem;
problemData.darcyData.p = problemData.p;
problemData.darcyData.qOrd = problemData.qOrd;

problemData.darcyData.t0 = problemData.t0;
problemData.darcyData.tEnd = problemData.tEnd;
problemData.darcyData.numSteps = problemData.numSteps;

if problemData.isCouplingDarcy
  problemData.darcyData.idCoupling = 3;
end % if

problemData.darcyData.isVisGrid = problemData.isVisGrid;
problemData.darcyData.isVisSol = problemData.isVisSol;
problemData.darcyData.outputFrequency = problemData.outputFrequency;
problemData.darcyData.outputBasename = ['output' filesep 'darcyVert_' problemData.outputBasename];
problemData.darcyData.outputTypes = problemData.outputTypes;

problemData.darcyData = problemData.darcySteps.configureProblem(problemData.darcyData);

%% Free-flow problem
problemData.sweData = struct;
problemData.sweData.problemName = 'sweVert';
problemData.sweData.testcase = problemData.testcase;

problemData.sweData.numElem = problemData.numElem;
problemData.sweData.p = problemData.p;
problemData.sweData.qOrd = problemData.qOrd;

problemData.sweData.t0 = problemData.t0;
problemData.sweData.tEnd = problemData.tEnd;
problemData.sweData.numSteps = problemData.numSteps * problemData.numSubSteps;

if problemData.isCouplingSWE
  problemData.sweData.idCoupling = 1;
end % if

problemData.sweData.isVisGrid = problemData.isVisGrid;
problemData.sweData.isVisSol = problemData.isVisSol;
problemData.sweData.outputFrequency = problemData.outputFrequency * problemData.numSubSteps;
problemData.sweData.outputBasename = ['output' filesep 'sweVert_' problemData.outputBasename];
problemData.sweData.outputTypes = problemData.outputTypes;

problemData.sweData = problemData.sweSteps.configureProblem(problemData.sweData);

%% Extract derived data
problemData.generateGrid1D = problemData.sweData.generateGrid1D;
end % function