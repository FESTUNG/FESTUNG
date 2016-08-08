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
function pd = configureProblem(pd)
%% Name of the problem
% Influences name of output files and specifies name of ADCIRC input files
pd = setdefault(pd, 'name', 'debug'); 

%% Configuration to use: 
% - 'debug' calls configureDebug()
% - 'analytical' calls configureAnalyticalTest()
% - 'ADCIRC' reads 'swe/fort_<name>.15'
pd = setdefault(pd, 'configSource', 'debug');

%% What kind of grid to use:
% - 'square' creates a unit square [0,1]x[0,1] with given pd.hmax,
%   open sea boundary in the east (type 4), and land boundary (type 1) on 
%   all other edges 
% - 'hierarchical' creates a unit square [0,1]x[0,1] with specified hmax
%   and performs uniform refinement according to parameter 'refinement'.
%   Boundary type 4 on east-boundary, 1 on all others.
% - 'ADCIRC' reads grid information from 'swe/fort_<name>.{14,17}'.
pd = setdefault(pd, 'gridSource', 'square');

%% Polynomial approximation order
% Piecewise constant (0), piecewise linear (1), or piecewise quadratic (2)
pd = setdefault(pd, 'p', 1);

%% Time stepping parameters
pd = setdefault(pd, 'schemeType', 'explicit'); % type of time stepping scheme ('explicit' or 'semi-implicit')

%% Model parameters
% Some may be overwritten by fort.15 config files
pd = setdefault(pd, 'typeFlux', 'Lax-Friedrichs'); % Type of interior flux ('Lax-Friedrichs', 'Roe')
pd = setdefault(pd, 'isRiemOS', true); % Riemann solver type on open sea boundary ('Lax-Friedrichs', 'Roe', or 'none')
pd = setdefault(pd, 'isRiemRiv', true); % Riemann solver type on river boundary ('Lax-Friedrichs', 'Roe', or 'none')
pd = setdefault(pd, 'typeBdrL', 'riemann'); % Flux type on land boundary ('reflected', 'natural', or 'riemann')
pd = setdefault(pd, 'averagingType', 'full-harmonic'); % Averaging type for variables when computing flux ('full-harmonic', 'semi-harmonic', 'mean')
pd = setdefault(pd, 'typeSlopeLim', 'linear'); % Slope limiter type ('linear', 'hierarch_vert', 'strict')
pd = setdefault(pd, 'slopeLimList', {}); % Apply slope limiter to specified variables ('h', 'uH', 'vH')
pd = setdefault(pd, 'isCoupling', false); % Compute velocity coefficients and flux of first unknown, e.g., for coupled transport problem
pd = setdefault(pd, 'elevTol', 20);

%% Visualization parameters
pd = setdefault(pd, 'isVisGrid', false); % Visualize computational grid
pd = setdefault(pd, 'isWaitbar', false); % Use waiting bar
pd = setdefault(pd, 'outputCount', 1); % Number of outputs over total simulation time
pd = setdefault(pd, 'outputTypes', {'vtk', 'tec'}); % Output file type
pd = setdefault(pd, 'outputList', { 'xi', 'h', 'z_b', 'f_c', 'velocity', 'momentum' }); % List of variables to visualize
pd = setdefault(pd, 'isVisStations', false); % Output stations

%% Simulation scenario specific parameters
switch pd.configSource
  case 'debug'
    pd = configureDebug(pd);
  case 'analytical'
    pd = configureAnalyticalTest(pd);
  case 'ADCIRC'
    pd = configureADCIRC(pd);
  otherwise
    error('Invalid config source.')
end % switch

end % function

%% Debugging
function pd = configureDebug(pd)
pd.isSolutionAvail = true;
pd.isRhsAvail = true;
pd.isTidalDomain = false;
pd.isHotstartInput = false;
pd.isHotstartOutput = false;
pd.schemeOrder = min(pd.p+1,3);

% Overwrite grid parameters
pd.gridSource = 'square';
pd.isSpherical = false;
pd.hmax = 2^-6; % Maximum element size of initial grid 

% Overwrite time-stepping parameters
pd.t0 = 0; % Start time of simulation
pd.numSteps = 300; % Number of time steps
pd.tEnd = pd.numSteps/3142*2*pi; % End time of simulation

pd.isAdaptiveTimestep = false; % Use adaptive timestep width
pd.dt = (pd.tEnd - pd.t0) / pd.numSteps;

pd.isSteadyState = false; % End simulation upon convergence

% Solution parameters
pd.gConst = 9.81;
pd.minTol = 0.001;

pd.isBottomFrictionNonlinear = true; % NOLIBF
pd.isBottomFrictionVarying = false; % NWP
pd.bottomFrictionCoef = 0;

% Ramping function, bathymetry, and Coriolis coefficient
pd.isRamp = false;
pd.ramp = @(t) 1;
pd.zbCont = @(x1,x2) -0.002*(x1==x1);
pd.fcCont = @(x1,x2) zeros(size(x1));

% Analytical solution
pd.xiCont = @(x1,x2,t) zeros(size(x1));
pd.uCont = @(x1,x2,t) 0.5 - x2;
pd.vCont = @(x1,x2,t) x1 - 0.5;

% Right hand side functions (derived from analytical solution)
pd.f0Cont = @(x1,x2,t) zeros(size(x1));
pd.f1Cont = @(x1,x2,t) 0.002*(0.5-x1);
pd.f2Cont = @(x1,x2,t) 0.002*(0.5-x2);

% Boundary conditions
pd.xiOSCont = @(x1,x2,t) pd.xiCont(x1,x2,t);
pd.isRivCont = true;
pd.xiRivCont = @(x1,x2,t) pd.xiCont(x1,x2,t);
pd.uRivCont = @(x1,x2,t) pd.uCont(x1,x2,t);
pd.vRivCont = @(x1,x2,t) pd.vCont(x1,x2,t);

% Hot-start output
pd.isHotStartOutput = false;
end % function

%% Analytical solution
function pd = configureAnalyticalTest(pd)
pd.isSolutionAvail = true;
pd.isRhsAvail = true;
pd.isTidalDomain = false;
pd.isHotstartInput = false;
pd.isHotstartOutput = false;
pd.schemeOrder = min(pd.p+1,3);

% Overwrite grid parameters
pd.gridSource = 'hierarchical';
pd.isSpherical = false; 
pd.hmax = 0.3; % Maximum element size of initial grid
pd.refinement = 6;  % Grid refinement level

% Overwrite time-stepping parameters
pd.t0 = 0; % Start time of simulation
pd.tEnd = 1; % End time of simulation
pd.numSteps = 400; % Number of time steps

pd.isAdaptiveTimestep = false; % Use adaptive timestep width
pd.dt = (pd.tEnd - pd.t0) / pd.numSteps;

pd.isSteadyState = false; % End simulation upon convergence

% Solution parameters
height = 0; % value of each component of bathymatry gradient
A = 0.01; % TODO height
B = 0.01;
C = 0.01;
pd.gConst = 9.81;
pd.minTol = 0.001;

pd.isBottomFrictionNonlinear = true; % NOLIBF
pd.isBottomFrictionVarying = false; % NWP
pd.bottomFrictionCoef = 0;

% Ramping function, bathymetry, and Coriolis coefficient
pd.isRamp = false;
pd.ramp = @(t) 1;
pd.zbCont = @(x1,x2) -0.1*(x1==x1); % TODO
pd.fcCont = @(x1,x2) 0*x1;

% % Analytical solution
% pd.xiCont = @(x1,x2,t) 0*x1;
% pd.uCont = @(x1,x2,t) x1.*(x1-1);
% pd.vCont = @(x1,x2,t) 0*x1;
% 
% % Auxiliary functions (derivatives etc.)
% pd.hCont = @(x1,x2,t) pd.xiCont(x1,x2,t) - pd.zbCont(x1,x2);
% pd.u_tCont = @(x1,x2,t) 0*x1;
% pd.u_xCont = @(x1,x2,t) 2*x1-1;
% pd.u_yCont = @(x1,x2,t) 0*x1;
% pd.v_tCont = @(x1,x2,t) 0*x1;
% pd.v_xCont = @(x1,x2,t) 0*x1;
% pd.v_yCont = @(x1,x2,t) 0*x1;
% pd.h_tCont = @(x1,x2,t) 0*x1;
% pd.h_xCont = @(x1,x2,t) -height*(x1==x1);
% pd.h_yCont = @(x1,x2,t) -height*(x1==x1);

% Analytical solution
pd.xiCont = @(x1,x2,t) 0*x1;
pd.uCont = @(x1,x2,t) A*sin(pi*x1);
pd.vCont = @(x1,x2,t) 0*x1;
% 
% Auxiliary functions (derivatives etc.)
pd.hCont = @(x1,x2,t) pd.xiCont(x1,x2,t) - pd.zbCont(x1,x2);
% pd.zb_xCont = @(x1,x2) height*(x1==x1); % TODO
% pd.zb_yCont = @(x1,x2) height*(x1==x1);
pd.u_tCont = @(x1,x2,t)  0*x1;
pd.u_xCont = @(x1,x2,t)  A*pi*cos(pi*x1);
pd.u_yCont = @(x1,x2,t)  0*x1;
pd.v_tCont = @(x1,x2,t)  0*x1;
pd.v_xCont = @(x1,x2,t)  0*x1;
pd.v_yCont = @(x1,x2,t)  0*x1;
pd.h_tCont = @(x1,x2,t)  0*x1;
pd.h_xCont = @(x1,x2,t)  0*x1;
pd.h_yCont = @(x1,x2,t)  0*x1;

% % Analytical solution
% pd.xiCont = @(x1,x2,t) C*(sin(0.5*pi*(x1-t)) + sin(0.5*pi*(x2-t)));
% pd.uCont = @(x1,x2,t) A*sin(0.5*pi*(x2-t)).*sin(pi*x1);
% pd.vCont = @(x1,x2,t) B*sin(0.5*pi*(x1-t)).*sin(pi*x2);
% % 
% % Auxiliary functions (derivatives etc.)
% pd.hCont = @(x1,x2,t) pd.xiCont(x1,x2,t) - pd.zbCont(x1,x2);
% % pd.zb_xCont = @(x1,x2) height*(x1==x1); % TODO
% % pd.zb_yCont = @(x1,x2) height*(x1==x1);
% pd.u_tCont = @(x1,x2,t) -0.5*pi*A*cos(0.5*pi*(x2-t)).*sin(pi*x1);
% pd.u_xCont = @(x1,x2,t)      pi*A*sin(0.5*pi*(x2-t)).*cos(pi*x1);
% pd.u_yCont = @(x1,x2,t)  0.5*pi*A*cos(0.5*pi*(x2-t)).*sin(pi*x1);
% pd.v_tCont = @(x1,x2,t) -0.5*pi*B*sin(0.5*pi*(x1-t)).*sin(pi*x2);
% pd.v_xCont = @(x1,x2,t)  0.5*pi*B*cos(0.5*pi*(x1-t)).*sin(pi*x2);
% pd.v_yCont = @(x1,x2,t)      pi*B*sin(0.5*pi*(x1-t)).*sin(pi*x2);
% pd.h_tCont = @(x1,x2,t) -0.5*pi*C*(cos(0.5*pi*(x1-t)) + cos(0.5*pi*(x2-t)));
% pd.h_xCont = @(x1,x2,t)  0.5*pi*C*cos(0.5*pi*(x1-t)) - height;
% pd.h_yCont = @(x1,x2,t)  0.5*pi*C*cos(0.5*pi*(x2-t)) - height;

% Right hand side functions (derived from analytical solution)
% This can be used for any solution, the user just has to specify all
% necessary auxiliary functions
pd.f0Cont = @(x1,x2,t) pd.h_tCont(x1,x2,t) + ...
                      (pd.u_xCont(x1,x2,t) + pd.v_yCont(x1,x2,t)) .* pd.hCont(x1,x2,t) + ...
                      pd.uCont(x1,x2,t) .* pd.h_xCont(x1,x2,t) + ...
                      pd.vCont(x1,x2,t) .* pd.h_yCont(x1,x2,t);
pd.f1Cont = @(x1,x2,t) pd.u_tCont(x1,x2,t) .* pd.hCont(x1,x2,t) + pd.uCont(x1,x2,t) .* pd.h_tCont(x1,x2,t) + ...
                       ( 2 * pd.uCont(x1,x2,t) .* pd.u_xCont(x1,x2,t) + pd.gConst * pd.h_xCont(x1,x2,t) + pd.u_yCont(x1,x2,t) .* pd.vCont(x1,x2,t) + pd.uCont(x1,x2,t) .* pd.v_yCont(x1,x2,t) ) .* pd.hCont(x1,x2,t) + ...
                       pd.uCont(x1,x2,t) .* pd.uCont(x1,x2,t) .* pd.h_xCont(x1,x2,t) + pd.uCont(x1,x2,t) .* pd.vCont(x1,x2,t) .* pd.h_yCont(x1,x2,t) + ...
                       pd.gConst * height * pd.hCont(x1,x2,t) + ...
                       pd.bottomFrictionCoef * sqrt( pd.uCont(x1,x2,t) .* pd.uCont(x1,x2,t) + pd.vCont(x1,x2,t) .* pd.vCont(x1,x2,t) ) .* pd.uCont(x1,x2,t) - ...
                       pd.fcCont(x1,x2) .* pd.vCont(x1,x2,t) .* pd.hCont(x1,x2,t);
pd.f2Cont = @(x1,x2,t) pd.v_tCont(x1,x2,t) .* pd.hCont(x1,x2,t) + pd.vCont(x1,x2,t) .* pd.h_tCont(x1,x2,t) + ...
                       ( pd.u_xCont(x1,x2,t) .* pd.vCont(x1,x2,t) + pd.uCont(x1,x2,t) .* pd.v_xCont(x1,x2,t) + 2 * pd.vCont(x1,x2,t) .* pd.v_yCont(x1,x2,t) + pd.gConst * pd.h_yCont(x1,x2,t) ) .* pd.hCont(x1,x2,t) + ...
                       pd.uCont(x1,x2,t) .* pd.vCont(x1,x2,t) .* pd.h_xCont(x1,x2,t) + pd.vCont(x1,x2,t) .* pd.vCont(x1,x2,t) .* pd.h_yCont(x1,x2,t) + ...
                       pd.gConst * height * pd.hCont(x1,x2,t) + ...
                       pd.bottomFrictionCoef * sqrt( pd.uCont(x1,x2,t) .* pd.uCont(x1,x2,t) + pd.vCont(x1,x2,t) .* pd.vCont(x1,x2,t) ) .* pd.vCont(x1,x2,t) + ...
                       pd.fcCont(x1,x2) .* pd.uCont(x1,x2,t) .* pd.hCont(x1,x2,t);
                     
% Boundary conditions
pd.xiOSCont = pd.xiCont;
pd.isRivCont = true;
pd.xiRivCont = @(x1,x2,t) pd.xiCont(x1,x2,t);
pd.uRivCont = @(x1,x2,t) pd.uCont(x1,x2,t);
pd.vRivCont = @(x1,x2,t) pd.vCont(x1,x2,t);

% Hot-start output
pd.isHotStartOutput = false;
end % function

%% ADCIRC
function pd = configureADCIRC(pd)
pd.isSolutionAvail = false;
pd.isRhsAvail = false;

% Verify input files exist
assert(exist(['swe/fort_' pd.name '.14'], 'file') == 2, ['Mesh file "swe/fort_' pd.name '.14" not found!'])
assert(exist(['swe/fort_' pd.name '.17'], 'file') == 2, ['Mesh file "swe/fort_' pd.name '.17" not found!'])
assert(exist(['swe/fort_' pd.name '.15'], 'file') == 2, ['Config file "swe/fort_' pd.name '.15" not found!'])

%% Read parameter file
h = getFunctionHandle('swe/readConfigADCIRC');
pd.configADCIRC = h(['swe/fort_' pd.name '.15']); % TODO getFunctionHandle oder execin?

%% Map ADCIRC variables to internal names
% Constants
pd.schemeOrder = pd.configADCIRC.IRK+1;
pd.minTol = pd.configADCIRC.H0;
pd.gConst = pd.configADCIRC.G;

% Simulation time
pd.t0 = pd.configADCIRC.STATIM;
pd.tEnd = pd.t0 + pd.configADCIRC.RNDAY * 86400;
pd.dt = pd.configADCIRC.DT;
pd.numSteps = round((pd.tEnd - pd.t0) / pd.dt);

% Adaptive time stepping
pd.isAdaptiveTimestep = pd.configADCIRC.NDTVAR == 1;

% Steady state simulation
pd.isSteadyState = pd.configADCIRC.ITRANS == 1;
pd.convergenceCriterion = pd.configADCIRC.CONVCR;

% Coordinate system
pd.isSpherical = pd.configADCIRC.ICS == 2;

% Bottom friction
pd.isBottomFrictionVarying = pd.configADCIRC.NWP == 1;
assert(~pd.isBottomFrictionVarying, 'Spatially varying bottom friction not implemented.');

pd.isBottomFrictionNonlinear = pd.configADCIRC.NOLIBF == 1;
if pd.isBottomFrictionNonlinear
  pd.bottomFrictionCoef = pd.configADCIRC.CF;
else
  pd.bottomFrictionCoef = pd.configADCIRC.TAU;
end % if

% Ramping
switch pd.configADCIRC.NRAMP
  case 0
    pd.isRamp = false;
    pd.ramp = @(t_days) 1;
  case 1
    pd.isRamp = true;
    pd.ramp = @(t_days) ((t_days - pd.configADCIRC.STATIM) / pd.configADCIRC.DRAMP) * ...
                        (t_days < pd.configADCIRC.STATIM + pd.configADCIRC.DRAMP) + ...
                        (t_days >= pd.configADCIRC.STATIM + pd.configADCIRC.DRAMP);
  otherwise
    error('Invalid value for NRAMP')
end % switch
    
% Newtonian tidal potential
pd.isTidalDomain = pd.configADCIRC.NTIP == 1;

% Boundary conditions
pd.isRivCont = false;

% Hot-start input and output
pd.isHotstartInput = pd.configADCIRC.IHOT == 1;
pd.hotstartInput = pd.configADCIRC.OUTP;

pd.isHotstartOutput = pd.configADCIRC.NHSTAR == 1;
pd.hotstartOutputFrequency = pd.configADCIRC.NHSINC;
end % function
