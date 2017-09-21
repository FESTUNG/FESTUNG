% Fills the problemData-struct with all basic configuration options.
% Problem parameters are to be modified inside this routine.

%===============================================================================
%> @file template/configureProblem.m
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
%> The struct must provide a positive numeric parameter <code>numSteps</code>
%> that specifies the number of iterative steps (e.g., number of time
%> steps) which are to be executed in the main part of the solution
%> algorithm.
%> For stationary problems set this to 1.
%>
%> @param  problemData  A (probably) empty struct with problem parameters.
%>                      @f$[\text{struct}]@f$
%>
%> @retval problemData  A struct with all necessary parameters and definitions
%>                      for the problem description. @f$[\text{struct}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Hennes Hajduk, Balthasar Reuter, Florian Frank, Vadym Aizinger
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

%% Configuration to use: 
% - 'rotation' calls configureRotation()
% - 'analytical' calls configureAnalyticalTest()
% - 'ADCIRC' reads 'swe/fort_<name>.15'
problemData = setdefault(problemData, 'configSource', 'ADCIRC');

%% What kind of grid to use:
% - 'square' creates a unit square [0,1]x[0,1] with given pd.hmax,
%   open sea boundary in the east (type 4), and land boundary (type 1) on 
%   all other edges 
% - 'hierarchical' creates a unit square [0,1]x[0,1] with specified hmax
%   and performs uniform refinement according to parameter 'refinement'.
%   Boundary type 4 on east-boundary, 1 on all others.
% - 'ADCIRC' reads grid information from 'swe/fort_<name>.{14,17}'.
problemData = setdefault(problemData, 'gridSource', 'ADCIRC');
problemData = setdefault(problemData, 'refinement', 0);

% Polynomial approximation order
problemData = setdefault(problemData, 'p', 1);

% Runge-Kutta order
problemData = setdefault(problemData, 'ordRK', min(problemData.p+1,3)); % as of now both models have to use the same RK method

% Time stepping specification, as of now both models have to use the same number of time steps
switch problemData.configSource
  case 'rotation'
    problemData = setdefault(problemData, 'hmax', 2^-6);
    problemData = setdefault(problemData, 'numSteps', 3124);
    problemData = setdefault(problemData, 'tEnd', (problemData.numSteps/3142)*2*pi);
  case 'biological'
    problemData.isSolutionAvailable = false;
    problemData = setdefault(problemData, 'hmax'      , 2^-6);  % maximum edge length of triangle
    problemData = setdefault(problemData, 'numSteps'  , 3142);  % number of time steps
    problemData = setdefault(problemData, 'tEnd'      , (problemData.numSteps/3142)*2*pi);  % end time
  case 'analytical'
    problemData = setdefault(problemData, 'hmax', 200);
    problemData = setdefault(problemData, 'numSteps', 200*2^(problemData.refinement+problemData.p));
    problemData = setdefault(problemData, 'tEnd', 500);
  case 'ADCIRC'
    problemData = setdefault(problemData, 'hmax', 1);
    problemData = setdefault(problemData, 'numSteps', 1555200);
    problemData = setdefault(problemData, 'tEnd', 7776000);
  otherwise
    error('Invalid config source.')
end % switch

% Configuration for shallow water solver
problemData.sweData = struct;

% Specification of configuration type for shallow water model
switch problemData.configSource
  case 'rotation'
    problemData.sweData.configSource = 'debug';
  case 'biological'
    problemData.sweData.configSource = 'debug';
  case 'analytical'
    problemData.sweData.configSource = 'analytical';
  case 'ADCIRC'
    problemData.isSolutionAvailable = false;
    problemData.name = 'galv';
    problemData.sweData.configSource = 'ADCIRC';
    problemData.sweData.name = problemData.name;
  otherwise
    error('Invalid config source.')
end % switch

problemData.sweData.isCoupling = true;

problemData.sweData.gridSource = problemData.gridSource;
problemData.sweData.refinement = problemData.refinement;
problemData.sweData.hmax = problemData.hmax;
problemData.sweData.p = problemData.p;
problemData.sweData.schemeOrder = problemData.ordRK;
problemData.sweData.tEnd = problemData.tEnd;
problemData.sweData.numSteps = problemData.numSteps;

h = getFunctionHandle('swe/configureProblem');
problemData.sweData = h(problemData.sweData);

problemData.ordRK = problemData.sweData.schemeOrder; % in case ADCIRC scheme order is used

% Configuration for transport solver
problemData.transportData = struct;

% Specification of configuration type for transport model and solution
switch problemData.configSource
  case 'rotation'
    problemData.transportData.configSource = 'rotation';
    problemData.transportData.hCont = @(t,x1,x2) problemData.sweData.xiCont(x1,x2,t) - problemData.sweData.zbCont(x1,x2);
  case 'biological'
    problemData.transportData.configSource = 'biological';
    problemData.transportData.hCont = @(t,x1,x2) problemData.sweData.xiCont(x1,x2,t) - problemData.sweData.zbCont(x1,x2);
  case 'analytical'
    problemData.transportData.configSource = 'analytical';
    
    % analytical functions from swe necessary for convergence test
    problemData.transportData.hCont = @(t,x1,x2) problemData.sweData.hCont(x1,x2,t);
    problemData.transportData.h_tCont = @(t,x1,x2) problemData.sweData.h_tCont(x1,x2,t);
    problemData.transportData.h_xCont = @(t,x1,x2) problemData.sweData.h_xCont(x1,x2,t);
    problemData.transportData.h_yCont = @(t,x1,x2) problemData.sweData.h_yCont(x1,x2,t);
    problemData.transportData.uCont = @(t,x1,x2) problemData.sweData.uCont(x1,x2,t);
    problemData.transportData.u_xCont = @(t,x1,x2) problemData.sweData.u_xCont(x1,x2,t);
    problemData.transportData.vCont = @(t,x1,x2) problemData.sweData.vCont(x1,x2,t);
    problemData.transportData.v_yCont = @(t,x1,x2) problemData.sweData.v_yCont(x1,x2,t);
  case 'ADCIRC'
    problemData.transportData.configSource = 'ADCIRC';
    problemData.transportData.name = problemData.name;
    problemData.transportData.domainADCIRC = getFunctionHandle('swe/domainADCIRC');
  otherwise
    error('Invalid config source.')
end % switch

problemData.transportData.isCoupling = true;

problemData.transportData.gridSource = problemData.gridSource;
problemData.transportData.refinement = problemData.refinement;
problemData.transportData.hmax = problemData.hmax;
problemData.transportData.p = problemData.p;
problemData.transportData.ordRK = problemData.ordRK;
problemData.transportData.tEnd = problemData.tEnd;
problemData.transportData.numSteps = problemData.numSteps;

problemData.transportData.isVisGrid = false; % visualization of grid

h = getFunctionHandle('transport/configureProblem');
problemData.transportData = h(problemData.transportData);
end % function