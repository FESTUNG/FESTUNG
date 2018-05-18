% Sets all problem parameters and fills the problemData-struct with all basic 
% configuration options.

%===============================================================================
%> @file
%>
%> @brief Sets all problem parameters and fills the problemData-struct with
%>        all basic configuration options.
%===============================================================================
%>
%> @brief Sets all problem parameters and fills the problemData-struct with
%>        all basic configuration options.
%>
%> This routine is called before any other function for the problem.
%> It defines all problem parameters and should be the only file users have
%> to adapt to their needs.
%>
%> Here, you can select the direction of the coupling, time stepping
%> parameterizations etc.
%> The actual problem definitions are in darcy_2dv/configureProblem.m or
%> swe_2dv/getTestcase.m
%>
%> For an outline of the models, see @link darcy_2dv @endlink and
%> @link swe_2dv @endlink.
%>
%> A detailed description can be found in @ref RRAFK2018.
%>
%> @param  problemData  A (probably) empty struct with problem parameters.
%>                      @f$[\text{struct}]@f$
%>
%> @retval problemData  A struct with all necessary parameters and definitions
%>                      for the problem description. @f$[\text{struct}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2018 Balthasar Reuter, Florian Frank, Vadym Aizinger
%>
%> @author Balthasar Reuter, 2018
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
problemData = setdefault(problemData, 'testcase', 'showcase');

% Enable coupling
problemData = setdefault(problemData, 'isCouplingDarcy', true);
problemData = setdefault(problemData, 'isCouplingSWE', true);

% Number of elements in x- and y-direction
problemData = setdefault(problemData, 'numElem', [32 8]);
problemData = setdefault(problemData, 'numElemDarcy', problemData.numElem);
problemData = setdefault(problemData, 'numElemSWE', problemData.numElem);

% Local polynomial approximation order (0 to 5)
problemData = setdefault(problemData, 'p', 1);
problemData = setdefault(problemData, 'pDarcy', problemData.p);
problemData = setdefault(problemData, 'pSWE', problemData.p);

% Order of quadrature rule
problemData = setdefault(problemData, 'qOrd', 2 * max(problemData.pDarcy, problemData.pSWE) + 1);

% Time stepping parameters
problemData = setdefault(problemData, 't0', 0);  % start time
problemData = setdefault(problemData, 'tEnd', 600);  % end time
problemData = setdefault(problemData, 'numSteps', ceil(problemData.tEnd/0.2));  % number of time steps
problemData = setdefault(problemData, 'numSubSteps', 10); % number of free-flow steps per sub-surface step

% Visualization settings
problemData = setdefault(problemData, 'isVisGrid', false);  % visualization of grid
problemData = setdefault(problemData, 'isVisSol', true);  % visualization of solution
problemData = setdefault(problemData, 'outputFrequency', 100); % no visualization of every timestep

couplingString = '_';
if problemData.isCouplingDarcy, couplingString = [couplingString 'couplingDarcy']; end
if problemData.isCouplingSWE, couplingString = [couplingString 'couplingSWE']; end

problemData = setdefault(problemData, 'outputBasename', [problemData.testcase couplingString]); 
problemData = setdefault(problemData, 'outputTypes', { 'vtk' });  % Type of visualization files ('vtk, 'tec')

%% Function handles for steps of the sub-problems
problemData.darcySteps = getStepHandles('darcy_2dv');
problemData.sweSteps = getStepHandles('swe_2dv');

%% Sub-surface problem
problemData.darcyData = struct;
problemData.darcyData.problemName = 'darcy_2dv';
problemData.darcyData.testcase = problemData.testcase;

problemData.darcyData.numElem = problemData.numElemDarcy;
problemData.darcyData.p = problemData.pDarcy;
problemData.darcyData.qOrdMax = problemData.qOrd;

problemData.darcyData.t0 = problemData.t0;
problemData.darcyData.tEnd = problemData.tEnd;
problemData.darcyData.numSteps = problemData.numSteps;

problemData.darcyData.isCoupling = problemData.isCouplingDarcy;
problemData.darcyData.isStationary = false;
problemData.darcyData.isJumpCoupling = true;

problemData.darcyData.isVisGrid = problemData.isVisGrid;
problemData.darcyData.isVisSol = problemData.isVisSol;
problemData.darcyData.outputFrequency = problemData.outputFrequency;
problemData.darcyData.outputBasename = ['output' filesep 'darcy_2dv_' problemData.outputBasename];
problemData.darcyData.outputTypes = problemData.outputTypes;

problemData.darcyData = problemData.darcySteps.configureProblem(problemData.darcyData);

%% Free-flow problem
problemData.sweData = struct;
problemData.sweData.problemName = 'swe_2dv';
problemData.sweData.testcase = problemData.testcase;

problemData.sweData.numElem = problemData.numElemSWE;
problemData.sweData.p = problemData.pSWE;
problemData.sweData.qOrdMax = problemData.qOrd;

problemData.sweData.t0 = problemData.t0;
problemData.sweData.tEnd = problemData.tEnd;
problemData.sweData.numSteps = problemData.numSteps * problemData.numSubSteps;

problemData.sweData.isCoupling = problemData.isCouplingSWE;

problemData.sweData.isVisGrid = problemData.isVisGrid;
problemData.sweData.isVisSol = problemData.isVisSol;
problemData.sweData.outputFrequency = problemData.outputFrequency * problemData.numSubSteps;
problemData.sweData.outputBasename = ['output' filesep 'swe_2dv_' problemData.outputBasename];
problemData.sweData.outputTypes = problemData.outputTypes;

problemData.sweData = problemData.sweSteps.configureProblem(problemData.sweData);

%% Extract derived data
problemData.generateGrid1D = problemData.sweData.generateGrid1D;
end % function
