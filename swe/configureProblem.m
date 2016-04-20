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
problemData.name = 'analytical'; % Name of the problem
problemData.p = 1; % Polynomial approximation order
problemData.hmax = 0.3; % Maximum element size of initial grid for analytical tests
problemData.refinement = 0;  % Grid refinement level

%% Visualization parameters
problemData.isVisGrid = false; % Visualize computational grid
problemData.isWaitbar = false; % Use waiting bar
problemData.outputCount = 20; % Number of outputs over total simulation time
problemData.outputTypes = 'vtk'; % Output file type
problemData.outputList = { 'u', 'uH', 'v', 'vH', 'xi', 'h', 'zb', 'fc' }; % List of variables to visualize

%% Time stepping parameters
problemData.scheme = 'explicit'; % type of time stepping scheme ('explicit' or 'implicit')
problemData.t0 = 0;
problemData.tEnd = 1;
problemData.numSteps = 150;

%% Model parameters
problemData.isOSRiem = true; % apply Riemann solver to open sea boundary
problemData.fluxType = 'Lax-Friedrichs';
problemData.averaging = 'full-harmonic';
problemData.minValueHeight = 1.e-3; % Minimum water level, values below are reset to this value

%% Analytical solution
problemData.isSolutionAvail = true;

% Solution parameters
height = 0.05;
A = 0.1;
B = 0.1;
C = 0.1;
problemData.gConst = 9.81;
problemData.bottomFric = 1;

% Ramping function, bathymetry, and Coriolis coefficient
problemData.ramp = @(t) 1;
problemData.zbCont = @(x1,x2) -height * (2 - x1 - x2);
problemData.fcCont = @(x1,x2) x1;

% Analytical solution
problemData.xiCont = @(x1,x2,t) C*(cos(0.5*pi*(x1-t)) + cos(0.5*pi*(x2-t))) - height*(2-x1-x2);
problemData.uCont = @(x1,x2,t) A*sin(pi*x1)*cos(2*pi*t);
problemData.vCont = @(x1,x2,t) B*sin(pi*x2)*cos(2*pi*t);

% Auxiliary functions (derivatives etc.)
problemData.hCont = @(x1,x2,t) problemData.xiCont(x1,x2,t) - problemData.zbCont(x1,x2);
problemData.u_tCont = @(x1,x2,t)   -2*A*pi*sin(pi*x1)*sin(2*pi*t);
problemData.u_xCont = @(x1,x2,t)      A*pi*cos(pi*x1)*cos(2*pi*t);
problemData.v_tCont = @(x1,x2,t)   -2*B*pi*sin(pi*x2)*sin(2*pi*t);
problemData.v_yCont = @(x1,x2,t)      B*pi*cos(pi*x2)*cos(2*pi*t);
problemData.h_tCont = @(x1,x2,t)  0.5*C*pi * ( sin(0.5*pi*(x1-t)) + sin(0.5*pi*(x2-t)) );
problemData.h_xCont = @(x1,x2,t) -0.5*C*pi *   sin(0.5*pi*(x1-t))                       ;
problemData.h_yCont = @(x1,x2,t) -0.5*C*pi *                       sin(0.5*pi*(x2-t))   ;

% Right hand side functions (derived from analytical solution)
problemData.f0Cont = @(x1,x2,t) problemData.h_tCont(x1,x2,t) + ...
                      (problemData.u_xCont(x1,x2,t) + problemData.v_yCont(x1,x2,t)) .* problemData.hCont(x1,x2,t) + ...
                      problemData.uCont(x1,x2,t) .* problemData.h_xCont(x1,x2,t) + ...
                      problemData.vCont(x1,x2,t) .* problemData.h_yCont(x1,x2,t);
problemData.f1Cont = @(x1,x2,t) problemData.u_tCont(x1,x2,t) .* problemData.hCont(x1,x2,t) + ...
                       problemData.uCont(x1,x2,t) .* problemData.h_tCont(x1,x2,t) + ...
                       ( 2 * problemData.uCont(x1,x2,t) .* problemData.u_xCont(x1,x2,t) + problemData.gConst * problemData.h_xCont(x1,x2,t) ) .* problemData.hCont(x1,x2,t) + ... 
                       problemData.uCont(x1,x2,t) .* problemData.uCont(x1,x2,t) .* problemData.h_xCont(x1,x2,t) + ...
                       problemData.uCont(x1,x2,t) .* problemData.v_yCont(x1,x2,t) .* problemData.hCont(x1,x2,t) + ... 
                       problemData.uCont(x1,x2,t) .* problemData.vCont(x1,x2,t) .* problemData.h_yCont(x1,x2,t) + ...
                       problemData.gConst * height * problemData.hCont(x1,x2,t) + ...
                       problemData.bottomFric * sqrt( problemData.uCont(x1,x2,t) .* problemData.uCont(x1,x2,t) + problemData.vCont(x1,x2,t) .* problemData.vCont(x1,x2,t) ) .* problemData.uCont(x1,x2,t) - ...
                       problemData.fcCont(x1,x2) .* problemData.vCont(x1,x2,t) .* problemData.hCont(x1,x2,t);
problemData.f2Cont = @(x1,x2,t) problemData.v_tCont(x1,x2,t) .* problemData.hCont(x1,x2,t) + ...
                       problemData.vCont(x1,x2,t) .* problemData.h_tCont(x1,x2,t) + ...
                       problemData.u_xCont(x1,x2,t) .* problemData.vCont(x1,x2,t) .* problemData.hCont(x1,x2,t) + ...
                       problemData.uCont(x1,x2,t) .* problemData.vCont(x1,x2,t) .* problemData.h_xCont(x1,x2,t) + ...
                       2 * problemData.vCont(x1,x2,t) .* problemData.v_yCont(x1,x2,t) .* problemData.hCont(x1,x2,t) + ...
                       problemData.vCont(x1,x2,t) .* problemData.vCont(x1,x2,t) .* problemData.h_yCont(x1,x2,t) + ... 
                       problemData.gConst * problemData.h_yCont(x1,x2,t) .* problemData.hCont(x1,x2,t) + ...
                       problemData.gConst * height * problemData.hCont(x1,x2,t) + ...
                       problemData.bottomFric * sqrt( problemData.uCont(x1,x2,t) .* problemData.uCont(x1,x2,t) + problemData.vCont(x1,x2,t) .* problemData.vCont(x1,x2,t) ) .* problemData.vCont(x1,x2,t) + ...
                       problemData.fcCont(x1,x2) .* problemData.uCont(x1,x2,t) .* problemData.hCont(x1,x2,t);
end % function
