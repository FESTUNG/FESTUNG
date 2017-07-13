% Fills the problemData-struct with all basic configuration options.
% Problem parameters are to be modified inside this routine.

%===============================================================================
%> @file advection/configureProblem.m
%>
%> @brief Fills the problemData-struct with all basic configuration options.
%>        Problem parameters are to be modified inside this routine.
%===============================================================================
%>
%> @brief Fills the problemData-struct with all basic configuration options.
%>        Problem parameters are to be modified inside this routine.
%>
%> This routine is called before any other function for the problem.
%> It defines all problem parameters and should be the only file users have
%> to adapt to their needs.
%> 
%> It provides all configuration options for the numerical solver to 
%> approximate solutions 
%> @f$c:\overline{J}\times\overline{\Omega}\rightarrow\mathbb{R}@f$ 
%> of the advection equation
%> @f{align*}{
%> \partial_t c  + \nabla\cdot (\mathbf{u}\,c) &\;=\; f            &&\text{in}~J\times\Omega\,,\\
%> c                                           &\;=\; c_\mathrm{D} &&\text{on}~J\times{\partial\Omega}_{\mathrm{in}}\,,\\
%> c                                           &\;=\; c^0          &&\text{on}~\{0\}\times\Omega\,.
%> @f}
%> The velocity @f$\mathbf{u}:J\times\Omega\rightarrow\mathbb{R}^2@f$ and 
%> right hand side@f$f:J\times\Omega\rightarrow \mathbb{R}@f$
%> may vary in time and space. 
%> A detailed description can be found in @ref RAWFK2016.
%>
%> Please read the inline-comments in the code for the meaning of each
%> configuration option.
%>
%> @param  problemData  A (probably) empty struct with problem parameters.
%>                      @f$[\text{struct}]@f$
%>
%> @retval problemData  A struct with all necessary parameters and definitions
%>                      for the problem description. @f$[\text{struct}]@f$
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
function problemData = configureProblem(problemData)
%% Parameters.
problemData = setdefault(problemData, 'testcase', 'LeVeque'); % 'SteadyProblem'); % 

problemData = setdefault(problemData, 'hmax', 2^-4); % maximum edge length of triangle
problemData = setdefault(problemData, 'p', 2); % local polynomial degree
problemData = setdefault(problemData, 'ordRK', min(problemData.p+1,4)); % order of Runge-Kutta time stepper.
problemData = setdefault(problemData, 'numSteps', 160); % number of time steps
problemData = setdefault(problemData, 'tEnd', 2*pi); % end time

problemData = setdefault(problemData, 'isVisGrid', false);
problemData = setdefault(problemData, 'isVisSol', false);

problemData = setdefault(problemData, 'outputFrequency', 20);
problemData = setdefault(problemData, 'outputBasename', ['output' filesep 'solution_hdg_advection']);
problemData = setdefault(problemData, 'outputTypes', {'vtk'});

%% HDG specific parameters
problemData = setdefault(problemData, 'stab', 1.0);
problemData = setdefault(problemData, 'isTrueLocalSolve', true);
problemData = setdefault(problemData, 'trueLocalSolveSize', 16); % 16 seems to be good in most cases

%% Testing?
problemData = setdefault(problemData, 'isInTesting', false);
problemData = setdefault(problemData, 'showWaitBar', false);
problemData = setdefault(problemData, 'showFprintfProgress', true);

%% Parameter check.
assert(problemData.p >= 0 && problemData.p <= 4, 'Polynomial order must be zero to four.')
assert(problemData.ordRK >= 1 && problemData.ordRK <= 4, 'Order of Runge Kutta must be zero to three.')
assert(problemData.hmax > 0, 'Maximum edge length must be positive.')
assert(problemData.numSteps > 0, 'Number of time steps must be positive.')

%% Coefficients and boundary data
problemData = execin([problemData.problemName filesep 'getTestcase'], problemData, problemData.testcase);


%% Domain and triangulation configuration.
% Triangulate unit square using pdetool (if available or Friedrichs-Keller otherwise).

% At the moment, the mesh will be defined by the test case

% if license('checkout','PDE_Toolbox')
%   problemData.generateGridData = @(hmax) domainPolygon([-0.5 0.5 0.5 -0.5], [-0.5 -0.5 0.5 0.5], hmax);
% else
%   fprintf('PDE_Toolbox not available. Using Friedrichs-Keller triangulation.\n');
%   problemData.generateGridData = @domainSquare;
% end % if
% Specify edge ids of boundary conditions

% Specify edge ids of boundary conditions
problemData.generateMarkE0Tint = @(g) g.idE0T == 0;
% Neumannboundaries are defined by test case
% problemData.generateMarkE0TbdrD = @(g) ~(g.markE0Tint | g.markE0TbdrN);

end % function
