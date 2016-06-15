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
pd.name = 'analytical_test'; 

%% Configuration to use: 
% - 'debug' calls configureDebug()
% - 'analytical' calls configureAnalyticalTest()
% - 'ADCIRC' reads 'swe/fort_<name>.15'
pd.configSource = 'analytical';

%% What kind of grid to use:
% - 'square' creates a unit square [0,1]x[0,1] with given pd.hmax,
%   open sea boundary in the east (type 4), and land boundary (type 1) on 
%   all other edges 
% - 'hierarchical' creates a unit square [0,1]x[0,1] with specified hmax
%   and performs uniform refinement according to parameter 'refinement'.
%   Boundary type 4 on east-boundary, 1 on all others.
% - 'ADCIRC' reads grid information from 'swe/fort_<name>.{14,17}'.
pd.gridSource = 'hierarchical';

%% Polynomial approximation order
% Piecewise constant (0), piecewise linear (1), or piecewise quadratic (2)
pd.p = 2;

%% Time stepping parameters
pd.schemeType = 'explicit'; % type of time stepping scheme ('explicit' or 'semi-implicit')
pd.schemeOrder = min(pd.p+1,2);

%% Model parameters
% Some may be overwritten by fort.15 config files
pd.typeFlux = 'Lax-Friedrichs'; % Type of interior flux ('Lax-Friedrichs', 'Roe')
pd.isRiemOS = true; % Riemann solver type on open sea boundary ('Lax-Friedrichs', 'Roe', or 'none')
pd.typeBdrL = 'riemann'; % Flux type on land boundary ('reflected', 'natural', or 'riemann')
pd.averagingType = 'full-harmonic'; % Averaging type for variables when computing flux ('full-harmonic', 'semi-harmonic', 'mean')
pd.typeSlopeLim = 'linear'; % Slope limiter type ('linear', 'hierarch_vert', 'strict')
pd.slopeLimList = {}; % Apply slope limiter to specified variables ('h', 'uH', 'vH')
pd.minTol = 0.001;

%% Visualization parameters
pd.isVisGrid = false; % Visualize computational grid
pd.isWaitbar = false; % Use waiting bar
pd.outputCount = 25; % Number of outputs over total simulation time
pd.outputTypes = 'vtk'; % Output file type
pd.outputList = { 'u', 'uH', 'v', 'vH', 'xi', 'h', 'zb', 'fc' }; % List of variables to visualize
pd.isVisStations = true; % Output stations

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
pd.isSolutionAvail = false;
pd.isRhsAvail = false;
pd.isTidalDomain = false;
pd.isSpherical = false;

% Overwrite grid parameters
pd.gridSource = 'square';
pd.isSpherical = false;
pd.hmax = 1; % Maximum element size of initial grid 

% Overwrite time-stepping parameters
pd.t0 = 0; % Start time of simulation
pd.tEnd = 0.01; % End time of simulation
pd.numSteps = 5; % Number of time steps
pd.dt = (pd.tEnd - pd.t0) / pd.numSteps;

% Solution parameters
pd.gConst = 9.81;

pd.isBottomFrictionNonlinear = true; % NOLIBF
pd.isBottomFrictionVarying = false; % NWP
pd.bottomFrictionCoef = 0;

% Ramping function, bathymetry, and Coriolis coefficient
pd.isRamp = false;
pd.ramp = @(t) 1;
pd.zbCont = @(x1,x2) -0.1*(1-x1<x2)-0.1;
pd.fcCont = @(x1,x2) zeros(size(x1));

% Analytical solution
pd.xiCont = @(x1,x2,t) zeros(size(x1));
pd.uCont = @(x1,x2,t) zeros(size(x1));
pd.vCont = @(x1,x2,t) zeros(size(x1));

% Boundary conditions
pd.xiOSCont = @(x1,x2,t) pd.xiCont(x1,x2,t);
pd.isRivCont = true;
pd.xiRivCont = @(x1,x2,t) pd.xiCont(x1,x2,t);
pd.uRivCont = @(x1,x2,t) pd.uCont(x1,x2,t);
pd.vRivCont = @(x1,x2,t) pd.vCont(x1,x2,t);
end % function

%% Analytical solution
function pd = configureAnalyticalTest(pd)
pd.isSolutionAvail = true;
pd.isRhsAvail = true;
pd.isTidalDomain = false;

% Overwrite grid parameters
pd.gridSource = 'hierarchical';
pd.isSpherical = false; 
pd.hmax = 0.3; % Maximum element size of initial grid
pd.refinement = 0;  % Grid refinement level

% Overwrite time-stepping parameters
pd.t0 = 0; % Start time of simulation
pd.tEnd = 1; % End time of simulation
pd.numSteps = 150; % Number of time steps
pd.dt = (pd.tEnd - pd.t0) / pd.numSteps;

% Solution parameters
height = 0.05;
A = 0.1;
B = 0.1;
C = 0.1;
pd.gConst = 9.81;

pd.isBottomFrictionNonlinear = true; % NOLIBF
pd.isBottomFrictionVarying = false; % NWP
pd.bottomFrictionCoef = 1;

% Ramping function, bathymetry, and Coriolis coefficient
pd.isRamp = false;
pd.ramp = @(t) 1;
pd.zbCont = @(x1,x2) -height * (2 - x1 - x2);
pd.fcCont = @(x1,x2) x1;

% Analytical solution
pd.xiCont = @(x1,x2,t) C*(cos(0.5*pi*(x1-t)) + cos(0.5*pi*(x2-t))) - height*(2-x1-x2);
pd.uCont = @(x1,x2,t) A*sin(pi*x1)*cos(2*pi*t);
pd.vCont = @(x1,x2,t) B*sin(pi*x2)*cos(2*pi*t);

% Auxiliary functions (derivatives etc.)
pd.hCont = @(x1,x2,t) pd.xiCont(x1,x2,t) - pd.zbCont(x1,x2);
pd.u_tCont = @(x1,x2,t)   -2*A*pi*sin(pi*x1)*sin(2*pi*t);
pd.u_xCont = @(x1,x2,t)      A*pi*cos(pi*x1)*cos(2*pi*t);
pd.v_tCont = @(x1,x2,t)   -2*B*pi*sin(pi*x2)*sin(2*pi*t);
pd.v_yCont = @(x1,x2,t)      B*pi*cos(pi*x2)*cos(2*pi*t);
pd.h_tCont = @(x1,x2,t)  0.5*C*pi * ( sin(0.5*pi*(x1-t)) + sin(0.5*pi*(x2-t)) );
pd.h_xCont = @(x1,x2,t) -0.5*C*pi *   sin(0.5*pi*(x1-t))                       ;
pd.h_yCont = @(x1,x2,t) -0.5*C*pi *                       sin(0.5*pi*(x2-t))   ;

% Right hand side functions (derived from analytical solution)
pd.f0Cont = @(x1,x2,t) pd.h_tCont(x1,x2,t) + ...
                      (pd.u_xCont(x1,x2,t) + pd.v_yCont(x1,x2,t)) .* pd.hCont(x1,x2,t) + ...
                      pd.uCont(x1,x2,t) .* pd.h_xCont(x1,x2,t) + ...
                      pd.vCont(x1,x2,t) .* pd.h_yCont(x1,x2,t);
pd.f1Cont = @(x1,x2,t) pd.u_tCont(x1,x2,t) .* pd.hCont(x1,x2,t) + ...
                       pd.uCont(x1,x2,t) .* pd.h_tCont(x1,x2,t) + ...
                       ( 2 * pd.uCont(x1,x2,t) .* pd.u_xCont(x1,x2,t) + pd.gConst * pd.h_xCont(x1,x2,t) ) .* pd.hCont(x1,x2,t) + ... 
                       pd.uCont(x1,x2,t) .* pd.uCont(x1,x2,t) .* pd.h_xCont(x1,x2,t) + ...
                       pd.uCont(x1,x2,t) .* pd.v_yCont(x1,x2,t) .* pd.hCont(x1,x2,t) + ... 
                       pd.uCont(x1,x2,t) .* pd.vCont(x1,x2,t) .* pd.h_yCont(x1,x2,t) + ...
                       pd.gConst * height * pd.hCont(x1,x2,t) + ...
                       pd.bottomFrictionCoef * sqrt( pd.uCont(x1,x2,t) .* pd.uCont(x1,x2,t) + pd.vCont(x1,x2,t) .* pd.vCont(x1,x2,t) ) .* pd.uCont(x1,x2,t) - ...
                       pd.fcCont(x1,x2) .* pd.vCont(x1,x2,t) .* pd.hCont(x1,x2,t);
pd.f2Cont = @(x1,x2,t) pd.v_tCont(x1,x2,t) .* pd.hCont(x1,x2,t) + ...
                       pd.vCont(x1,x2,t) .* pd.h_tCont(x1,x2,t) + ...
                       pd.u_xCont(x1,x2,t) .* pd.vCont(x1,x2,t) .* pd.hCont(x1,x2,t) + ...
                       pd.uCont(x1,x2,t) .* pd.vCont(x1,x2,t) .* pd.h_xCont(x1,x2,t) + ...
                       2 * pd.vCont(x1,x2,t) .* pd.v_yCont(x1,x2,t) .* pd.hCont(x1,x2,t) + ...
                       pd.vCont(x1,x2,t) .* pd.vCont(x1,x2,t) .* pd.h_yCont(x1,x2,t) + ... 
                       pd.gConst * pd.h_yCont(x1,x2,t) .* pd.hCont(x1,x2,t) + ...
                       pd.gConst * height * pd.hCont(x1,x2,t) + ...
                       pd.bottomFrictionCoef * sqrt( pd.uCont(x1,x2,t) .* pd.uCont(x1,x2,t) + pd.vCont(x1,x2,t) .* pd.vCont(x1,x2,t) ) .* pd.vCont(x1,x2,t) + ...
                       pd.fcCont(x1,x2) .* pd.uCont(x1,x2,t) .* pd.hCont(x1,x2,t);
                     
% Boundary conditions
pd.xiOSCont = pd.xiCont;
pd.isRivCont = false;
end % function

%% ADCIRC
function pd = configureADCIRC(pd)
pd.isSolutionAvail = false;
pd.isRhsAvail = false;
pd.isTidalDomain = false;

% Verify input files exist
assert(exist(['swe/fort_' pd.name '.14'], 'file') == 2, ['Mesh file "swe/fort_' pd.name '.14" not found!'])
assert(exist(['swe/fort_' pd.name '.17'], 'file') == 2, ['Mesh file "swe/fort_' pd.name '.17" not found!'])
assert(exist(['swe/fort_' pd.name '.15'], 'file') == 2, ['Config file "swe/fort_' pd.name '.15" not found!'])

%% Read parameter file
pd.configADCIRC = readConfigADCIRC(['swe/fort_' pd.name '.15']);

%% Map ADCIRC variables to internal names
% Constants
pd.gConst = pd.configADCIRC.G;

% Simulation time
pd.t0 = pd.configADCIRC.STATIM;
pd.tEnd = pd.configADCIRC.RNDAY * 86400;
pd.dt = pd.configADCIRC.DT;
pd.numSteps = round((pd.tEnd - pd.t0) / pd.dt);

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
end % function