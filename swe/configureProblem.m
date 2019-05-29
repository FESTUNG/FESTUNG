% Fills the problemData-struct with all basic configuration options.
% Problem parameters are to be modified inside this routine.

%===============================================================================
%> @file
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
pd = setdefault(pd, 'name', 'analytical');

%% Configuration to use: 
% - 'debug' calls configureDebug()
% - 'analytical' calls configureAnalyticalTest()
% - 'ADCIRC' reads 'swe/fort_<name>.15'
% - 'manual' calls configureManualADCIRC()
pd = setdefault(pd, 'configSource', 'analytical');

%% What kind of grid to use:
% - 'square' creates a unit square [0,1]x[0,1] with given pd.hmax,
%   open sea boundary in the east (type 4), and land boundary (type 1) on 
%   all other edges 
% - 'hierarchical' creates a unit square [0,100]x[0,100] with specified 
%   hmax and performs uniform refinement according to parameter 
%   'refinement'. All boundaries are type 3, i.e river.
% - 'ADCIRC' reads grid information from 'swe/fort_<name>.{14,17}'.
% - 'periodicRectangle' creates a structured grid over a rectangle with
%   periodic boundary conditions on all boundaries
% - 'load' reads grid information from the mat-file swe/grids/<name>.mat
pd = setdefault(pd, 'gridSource', 'hierarchical');

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
pd = setdefault(pd, 'averagingType', 'semi-harmonic'); % Averaging type for variables when computing flux ('full-harmonic', 'semi-harmonic', 'mean')
pd = setdefault(pd, 'slopeLimName', 'sequentialSegregated'); % Different approches for slope limiting: 'componentwise', 'sequentialSegregated', 'sequentialVector', 'vector', suffix '+norm'
pd = setdefault(pd, 'variant', 'transformation'); % Different variants used to obtain a rotation matrix for vector limiting, only in cases 'sequentialVector' or 'vector'
pd = setdefault(pd, 'angle', 0); % Angle used for backtransformation in vector slope limiting, only in cases 'sequentialVector' or 'vector' and if variant is 'fixed'
pd = setdefault(pd, 'typeSlopeLim', 'linear'); % Slope limiter type ('linear', 'hierarch_vert', 'strict')
pd = setdefault(pd, 'slopeLimList', {'height', 'velocity'}); % Apply slope limiter to specified variables ('height', 'momentum' 'velocity') 'momentum' and 'velocity' are exclusive
pd = setdefault(pd, 'isCoupling', false); % Compute velocity coefficients and flux of first unknown, e.g., for coupled transport problem
pd = setdefault(pd, 'elevTol', 20); % maximum absolute value for elevation

%% Visualization parameters
pd = setdefault(pd, 'isVisGrid', false); % Visualize computational grid
pd = setdefault(pd, 'isWaitbar', false); % Use waitbar
pd = setdefault(pd, 'outputTypes', cellstr('vtk')); % Output file type
pd = setdefault(pd, 'outputList', {'elevation', 'velocity'}); % List of variables to visualize
pd = setdefault(pd, 'isVisStations', false); % Output stations
pd = setdefault(pd, 'isOutputGrid', false); % Specify whether the grid should be saved as a mat-file

%% Simulation scenario specific parameters
switch pd.configSource
  case 'debug'
    pd = configureDebug(pd);
  case 'analytical'
    pd = configureAnalyticalTest(pd);
  case 'ADCIRC'
    pd = configureADCIRC(pd);
  case 'manual'
    pd = configureManualADCIRC(pd);
  otherwise
    error('Invalid config source.')
end % switch

end % function

%% Debugging
function pd = configureDebug(pd)
pd.isSolutionAvail = true;
pd.isRhsAvail = false;
pd.isTidalDomain = false;
pd.isHotstartInput = false;
pd.isHotstartOutput = true;
pd.hotstartOutputFrequency = 20000;
pd = setdefault(pd, 'schemeOrder', min(pd.p+1,3));
pd.isSlopeLim = true;

% Overwrite grid parameters
pd.X1 = [0 1000];
pd.X2 = [0 1];
pd.isSpherical = false;
pd.hx = 1;
pd.hy = 1;

% Overwrite time-stepping parameters
pd.t0 = 0; % Start time of simulation
pd.outputStart = 0 * ones(1,4);
pd = setdefault(pd, 'numSteps', 20000);  % number of time steps
pd = setdefault(pd, 'tEnd', 20);  % end time
pd = setdefault(pd, 'outputCount', 1); % Number of outputs over total simulation time

pd.isAdaptiveTimestep = false; % Use adaptive timestep width
pd.dt = (pd.tEnd - pd.t0) / pd.numSteps;

pd.isSteadyState = false; % End simulation upon convergence

% Solution parameters
pd.gConst = 9.81;
pd.minTol = 0.001;

pd.isBottomFrictionNonlinear = false; % NOLIBF
pd.isBottomFrictionVarying = false; % NWP
pd.bottomFrictionCoef = 0;

% Ramping function, bathymetry, and Coriolis coefficient
pd.isRamp = false;
pd.ramp = @(t) 1;
pd.zbCont = @(x1,x2) -(x1==x1);
pd.fcCont = @(x1,x2) 0*x1;

% Analytical solution
cm = 6.23416;
xA = @(t) 500 - t*sqrt(10*pd.gConst);
xB = @(t) 500 + t*(2*sqrt(10*pd.gConst) - 3*cm);
xC = @(t) 500 + t*(2*cm^2*(sqrt(10*pd.gConst) - cm))/(cm^2 - pd.gConst);
rot = @(x1,x2) x1; % sqrt(x1.^2+x2.^2);
pd.xiCont = @(x1,x2,t) 9*(rot(x1,x2)<xA(t)) + setNaN2Zero( (4/9/pd.gConst*(sqrt(10*pd.gConst) - (rot(x1,x2)-500)/2/t).^2 - 1) .* (rot(x1,x2)>=xA(t)) .* (rot(x1,x2)<xB(t)) ) ...
                        + (cm^2/pd.gConst - 1) .* (rot(x1,x2)>=xB(t)) .* (rot(x1,x2)<xC(t));
pd.uCont = @(x1,x2,t) setNaN2Zero( 2/3*((rot(x1,x2)-500)/t + sqrt(10*9.81)) .* (rot(x1,x2) >= xA(t)) .* (rot(x1,x2) < xB(t)) ...
                        + 2 * (sqrt(10*9.81) - cm)  .* (rot(x1,x2) >= xB(t)) .* (rot(x1,x2) < xC(t)) );
pd.vCont = @(x1,x2,t) 0*x1;

% Right hand side functions (derived from analytical solution)
pd.f0Cont = @(x1,x2,t) 0*x1;
pd.f1Cont = @(x1,x2,t) 0*x1;
pd.f2Cont = @(x1,x2,t) 0*x1;

% Boundary conditions
pd.isOSCont = true;
pd.xiOSCont = @(x1,x2,t) pd.xiCont(x1,x2,t);
pd.isRivCont = false;
pd.xiRivCont = @(x1,x2,t) pd.xiCont(x1,x2,t);
pd.uRivCont = @(x1,x2,t) pd.uCont(x1,x2,t);
pd.vRivCont = @(x1,x2,t) pd.vCont(x1,x2,t);

end % function

%% Analytical solution
function pd = configureAnalyticalTest(pd)
pd.isSolutionAvail = true;
pd.isRhsAvail = true;
pd.isTidalDomain = false;
pd.isHotstartInput = false;
pd.isHotstartOutput = false;
pd = setdefault(pd, 'schemeOrder', min(pd.p+1,3));
pd.isSlopeLim = true;

% Overwrite grid parameters
pd.gridSource = 'hierarchical';
pd.isSpherical = false; 
pd = setdefault(pd, 'refinement', 0);
pd = setdefault(pd, 'hmax', 10);

% Overwrite time-stepping parameters
pd.t0 = 0; % Start time of simulation
pd = setdefault(pd, 'numSteps', 100*2^(pd.refinement+pd.p));  % number of time steps
pd = setdefault(pd, 'tEnd', 100);  % end time
pd = setdefault(pd, 'outputCount', 10);  % Number of outputs over total simulation time

pd.isAdaptiveTimestep = false; % Use adaptive timestep width
pd.dt = (pd.tEnd - pd.t0) / pd.numSteps;

pd.isSteadyState = false; % End simulation upon convergence

pd.gConst = 1;
pd.minTol = 0.001;

pd.isBottomFrictionNonlinear = false; % NOLIBF
pd.isBottomFrictionVarying = false; % NWP
pd.bottomFrictionCoef = 0.0;

% Ramping function, bathymetry, and Coriolis coefficient
pd.isRamp = false;
pd.ramp = @(t) 1;
pd.zbCont = @(x1,x2) -(x1==x1);
pd.fcCont = @(x1,x2) 0*x1;

M = 1/2;
c1 = -0.04;
c2 = 0.02;
a = pi/6;
x0 = 30;
y0 = 40;
fCont = @(x,y,t) -c2 * ((x-x0 - M*t*cos(a)).^2 + (y-y0 - M*t*sin(a)).^2);
f_tCont = @(x,y,t) 2*M*c2 * ((x-x0 - M*t*cos(a))*cos(a) + (y-y0 - M*t*sin(a))*sin(a));
f_xCont = @(x,y,t) -2*c2 * (x-x0 - M*t*cos(a));
f_yCont = @(x,y,t) -2*c2 * (y-y0 - M*t*sin(a));

% Analytical solution: For this example the right hand sides are zero
pd.xiCont = @(x1,x2,t) -c1^2/(4*c2*pd.gConst) * exp(2*fCont(x1,x2,t));
pd.uCont = @(x1,x2,t) M*cos(a) + c1 * (x2-y0 - M*t*sin(a)) .* exp(fCont(x1,x2,t));
pd.vCont = @(x1,x2,t) M*sin(a) - c1 * (x1-x0 - M*t*cos(a)) .* exp(fCont(x1,x2,t));

% Auxiliary functions (derivatives etc.)
pd.hCont = @(x1,x2,t) pd.xiCont(x1,x2,t) - pd.zbCont(x1,x2);
pd.zb_xCont = @(x1,x2) 0*x1;
pd.zb_yCont = @(x1,x2) 0*x1;
pd.u_tCont = @(x1,x2,t)  c1 * (x2-y0 - M*t*sin(a)) .* exp(fCont(x1,x2,t)) .* f_tCont(x1,x2,t) - c1*M*sin(a) * exp(fCont(x1,x2,t));
pd.u_xCont = @(x1,x2,t)  c1 * (x2-y0 - M*t*sin(a)) .* exp(fCont(x1,x2,t)) .* f_xCont(x1,x2,t);
pd.u_yCont = @(x1,x2,t)  c1 * (x2-y0 - M*t*sin(a)) .* exp(fCont(x1,x2,t)) .* f_yCont(x1,x2,t) + c1 * exp(fCont(x1,x2,t));
pd.v_tCont = @(x1,x2,t) -c1 * (x1-x0 - M*t*cos(a)) .* exp(fCont(x1,x2,t)) .* f_tCont(x1,x2,t) + c1*M*cos(a) * exp(fCont(x1,x2,t));
pd.v_xCont = @(x1,x2,t) -c1 * (x1-x0 - M*t*cos(a)) .* exp(fCont(x1,x2,t)) .* f_xCont(x1,x2,t) - c1 * exp(fCont(x1,x2,t));
pd.v_yCont = @(x1,x2,t) -c1 * (x1-x0 - M*t*cos(a)) .* exp(fCont(x1,x2,t)) .* f_yCont(x1,x2,t);
pd.h_tCont = @(x1,x2,t) 2 * pd.xiCont(x1,x2,t) .* f_tCont(x1,x2,t);
pd.h_xCont = @(x1,x2,t) 2 * pd.xiCont(x1,x2,t) .* f_xCont(x1,x2,t);
pd.h_yCont = @(x1,x2,t) 2 * pd.xiCont(x1,x2,t) .* f_yCont(x1,x2,t);

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
                       pd.gConst * pd.zb_xCont(x1,x2) .* pd.hCont(x1,x2,t) + ...
                       (1-pd.isBottomFrictionNonlinear) * pd.bottomFrictionCoef * pd.uCont(x1,x2,t) .* pd.hCont(x1,x2,t) + ...
                       pd.isBottomFrictionNonlinear * pd.bottomFrictionCoef * sqrt( pd.uCont(x1,x2,t) .* pd.uCont(x1,x2,t) + pd.vCont(x1,x2,t) .* pd.vCont(x1,x2,t) ) .* pd.uCont(x1,x2,t) - ...
                       pd.fcCont(x1,x2) .* pd.vCont(x1,x2,t) .* pd.hCont(x1,x2,t);
pd.f2Cont = @(x1,x2,t) pd.v_tCont(x1,x2,t) .* pd.hCont(x1,x2,t) + pd.vCont(x1,x2,t) .* pd.h_tCont(x1,x2,t) + ...
                       ( pd.u_xCont(x1,x2,t) .* pd.vCont(x1,x2,t) + pd.uCont(x1,x2,t) .* pd.v_xCont(x1,x2,t) + 2 * pd.vCont(x1,x2,t) .* pd.v_yCont(x1,x2,t) + pd.gConst * pd.h_yCont(x1,x2,t) ) .* pd.hCont(x1,x2,t) + ...
                       pd.uCont(x1,x2,t) .* pd.vCont(x1,x2,t) .* pd.h_xCont(x1,x2,t) + pd.vCont(x1,x2,t) .* pd.vCont(x1,x2,t) .* pd.h_yCont(x1,x2,t) + ...
                       pd.gConst * pd.zb_yCont(x1,x2) .* pd.hCont(x1,x2,t) + ...
                       (1-pd.isBottomFrictionNonlinear) * pd.bottomFrictionCoef * pd.vCont(x1,x2,t) .* pd.hCont(x1,x2,t) + ...
                       pd.isBottomFrictionNonlinear * pd.bottomFrictionCoef * sqrt( pd.uCont(x1,x2,t) .* pd.uCont(x1,x2,t) + pd.vCont(x1,x2,t) .* pd.vCont(x1,x2,t) ) .* pd.vCont(x1,x2,t) + ...
                       pd.fcCont(x1,x2) .* pd.uCont(x1,x2,t) .* pd.hCont(x1,x2,t);
                     
% Boundary conditions
pd.isOSCont = false;
pd.xiOSCont = pd.xiCont;
pd.isRivCont = true;
pd.xiRivCont = pd.xiCont;
pd.uRivCont = pd.uCont;
pd.vRivCont = pd.vCont;

end % function

%% ADCIRC
function pd = configureADCIRC(pd)
pd.isSolutionAvail = false;
pd.isRhsAvail = false;

% Verify input files exist
assert(exist(['swe/fortFiles/' pd.name '.14'], 'file') == 2, ['Mesh file "swe/fortFiles/' pd.name '.14" not found!'])
assert(exist(['swe/fortFiles/' pd.name '.17'], 'file') == 2, ['Mesh file "swe/fortFiles/' pd.name '.17" not found!'])
assert(exist(['swe/fortFiles/' pd.name '.15'], 'file') == 2, ['Config file "swe/fortFiles/' pd.name '.15" not found!'])

%% Read parameter file
h = getFunctionHandle('swe/readConfigADCIRC');
pd.configADCIRC = h(['swe/fortFiles/' pd.name '.15']);

%% Map ADCIRC variables to internal names
% Constants
pd.p = pd.configADCIRC.IRK;
pd.schemeOrder = min(pd.configADCIRC.IRK+1, 3);

pd.isSlopeLim = mod(pd.configADCIRC.ISLOPE, 2) == 1;

pd.minTol = pd.configADCIRC.H0;
pd.gConst = pd.configADCIRC.G;

% Simulation time
pd.t0 = pd.configADCIRC.STATIM * 86400;
pd.tEnd = pd.t0 + pd.configADCIRC.RNDAY * 86400;
pd.dt = pd.configADCIRC.DT;
pd.numSteps = round((pd.tEnd - pd.t0) / pd.dt);
pd = setdefault(pd, 'outputCount', 180); % Number of outputs over total simulation time

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

pd.isVisStations = pd.configADCIRC.NOUTE || pd.configADCIRC.NOUTV;

% Boundary conditions
pd.isOSCont = false;
pd.isRivCont = false;

if pd.configADCIRC.NOUTGE
  pd.outputList = [pd.outputList, 'elevation'];
else
  pd.outputList = pd.outputList( setdiff( 1:length(pd.outputList), [find(strcmp(pd.outputList, 'elevation')), find(strcmp(pd.outputList, 'height'))] ) );
end % if
if pd.configADCIRC.NOUTGV
  pd.outputList = [pd.outputList, 'velocity'];
else
  pd.outputList = pd.outputList( setdiff( 1:length(pd.outputList), [find(strcmp(pd.outputList, 'velocity')), find(strcmp(pd.outputList, 'momentum'))] ) );
end % if
pd.outputList = unique(pd.outputList);

pd.outputStart = pd.t0 + 86400 * [ pd.configADCIRC.TOUTSGE, pd.configADCIRC.TOUTSE, pd.configADCIRC.TOUTSGV, pd.configADCIRC.TOUTSV ];
pd.outputEnd = pd.t0 + 86400 * [ pd.configADCIRC.TOUTFGE, pd.configADCIRC.TOUTFE, pd.configADCIRC.TOUTFGV, pd.configADCIRC.TOUTFV ];
pd.outputFrequency = [ pd.configADCIRC.NSPOOLGE, pd.configADCIRC.NSPOOLE, pd.configADCIRC.NSPOOLGV, pd.configADCIRC.NSPOOLV ];

% Hot-start input and output
pd.isHotstartInput = pd.configADCIRC.IHOT == 1;
pd.hotstartInput = pd.configADCIRC.OUTP;

pd.isHotstartOutput = pd.configADCIRC.NHSTAR == 1;
pd.hotstartOutputFrequency = pd.configADCIRC.NHSINC;
end % function

function pd = configureManualADCIRC(pd)

% Verify input files exist
assert(exist(['swe/fortFiles/' pd.name '.14'], 'file') == 2, ['Mesh file "swe/fortFiles/' pd.name '.14" not found!'])
assert(exist(['swe/fortFiles/' pd.name '.17'], 'file') == 2, ['Mesh file "swe/fortFiles/' pd.name '.17" not found!'])

% Overwrite grid config source
pd.gridSource = 'ADCIRC';

switch pd.name
  case 'bahamas'
    pd.isSolutionAvail = false;
    pd.isRhsAvail = false;
    pd.isTidalDomain = false;
    
    pd.isHotstartInput = false;
    pd.isHotstartOutput = false;
    pd = setdefault(pd, 'schemeOrder', min(pd.p+1,3));
    pd.isSlopeLim = false;

    % Overwrite grid parameters
    pd.isSpherical = false;
    pd.configADCIRC.SLAM0 = 0;
    pd.configADCIRC.SFEA0 = 0;

    % Overwrite time-stepping parameters
    pd = setdefault(pd, 't0', 0); % Start time of simulation
    pd = setdefault(pd, 'tEnd', 1036800);  % end time
    pd = setdefault(pd, 'dt', 15); % Time step size
    pd = setdefault(pd, 'numSteps', round((pd.tEnd - pd.t0) / pd.dt));
    pd = setdefault(pd, 'outputCount', 240); % Number of outputs over total simulation time

    pd.isAdaptiveTimestep = false; % Use adaptive timestep width

    pd.isSteadyState = false; % End simulation upon convergence

    % Solution parameters
    pd.gConst = 9.81;
    pd.minTol = 0.001;

    pd.isBottomFrictionNonlinear = true; % NOLIBF
    pd.isBottomFrictionVarying = false; % NWP
    pd.bottomFrictionCoef = 0.0090;
    
    % Ramping function
    pd.isRamp = false;
    pd.ramp = @(t_days) 1;

    % Coriolis coefficient
    pd.configADCIRC.NCOR = 0;
    pd.configADCIRC.CORI = 3.19e-5;

    % Boundary conditions
    pd.isOSCont = true;
    pd.xiOSCont = @(x1,x2,t)  ( 0.075 * cos(0.000067597751162*t - pi/180*194.806 ) ...
                              + 0.095 * cos(0.000072921165921*t - pi/180*206.265 ) ...
                              + 0.10  * cos(0.000137879713787*t - pi/180*340.0   ) ...
                              + 0.395 * cos(0.000140518917083*t - pi/180*  0.0   ) ...
                              + 0.06  * cos(0.000145444119418*t - pi/180*42.97180) ) * (x1==x1);

    pd.configADCIRC.NBFR = 0;

    pd.isRivCont = false;
  case 'test2'
    pd.isSolutionAvail = false;
    pd.isRhsAvail = false;
    pd.isTidalDomain = false;
    
    pd.isHotstartInput = true;
    pd.hotstartInput = 'swe/initialConditions/test2_lin.mat';
    pd.isHotstartOutput = false;
    
    pd = setdefault(pd, 'schemeOrder', min(pd.p+1,3));
    pd.isSlopeLim = true;
    
    % Overwrite grid parameters
    pd.isSpherical = false;
    pd.configADCIRC.SLAM0 = 0;
    pd.configADCIRC.SFEA0 = 0;
    
    % Overwrite time-stepping parameters
    pd = setdefault(pd, 't0', 0); % Start time of simulation
    pd = setdefault(pd, 'tEnd', 259.2);  % end time
    pd = setdefault(pd, 'dt', 0.1); % Time step size
    pd = setdefault(pd, 'numSteps', round((pd.tEnd - pd.t0) / pd.dt));
    pd = setdefault(pd, 'outputCount', 72); % Number of outputs over total simulation time

    pd.isAdaptiveTimestep = false; % Use adaptive timestep width
    
    % Steady state simulation
    pd.isSteadyState = true;
    pd.convergenceCriterion = 0.0031;
    
    % Solution parameters
    pd.gConst = 0.16;
    pd.minTol = 0.05;
    
    pd.isBottomFrictionNonlinear = false; % NOLIBF
    pd.isBottomFrictionVarying = false; % NWP
    pd.bottomFrictionCoef = 0.0;
    
    % Ramping function
    pd.isRamp = false;
    pd.ramp = @(t_days) 1;
    
    % Coriolis coefficient
    pd.configADCIRC.NCOR = 0;
    pd.configADCIRC.CORI = 0.0;
    
    pd.configADCIRC.NBFR = 0;
    pd.isOSCont = false;
    
    pd.isRivCont = true;
    pd.xiRivCont = @(x1,x2,t) 0*x1;
    pd.uRivCont = @(x1,x2,t) x1==x1;
    pd.vRivCont = @(x1,x2,t) 0*x1;
  case 'galv'
    pd.isSolutionAvail = false;
    pd.isRhsAvail = false;
    pd.isTidalDomain = false;
    
    pd.isHotstartInput = true;
    pd.hotstartInput = 'swe/initialConditions/galv10dStartup_lin.mat';
    pd.isHotstartOutput = false;
    pd = setdefault(pd, 'schemeOrder', min(pd.p+1,3));
    pd.isSlopeLim = false;

    % Overwrite grid parameters
    pd.isSpherical = false;
    pd.configADCIRC.SLAM0 = 0;
    pd.configADCIRC.SFEA0 = 0;

    % Overwrite time-stepping parameters
    pd = setdefault(pd, 't0', 864000); % Start time of simulation
    pd = setdefault(pd, 'tEnd', 1728000);  % end time
    pd = setdefault(pd, 'dt', 5); % Time step size
    pd = setdefault(pd, 'numSteps', round((pd.tEnd - pd.t0) / pd.dt));
    pd = setdefault(pd, 'outputCount', 200); % Number of outputs over total simulation time

    pd.isAdaptiveTimestep = false; % Use adaptive timestep width

    pd.isSteadyState = false; % End simulation upon convergence

    % Solution parameters
    pd.gConst = 9.81;
    pd.minTol = 0.25;

    pd.isBottomFrictionNonlinear = true; % NOLIBF
    pd.isBottomFrictionVarying = false; % NWP
    pd.bottomFrictionCoef = 0.0040;
    
    % Ramping function
    pd.isRamp = false;
    pd.ramp = @(t_days) 1;

    % Coriolis coefficient
    pd.configADCIRC.NCOR = 0;
    pd.configADCIRC.CORI = 7.07E-5;

    % Boundary conditions
    pd.isOSCont = true;
    pd.xiOSCont = @(x1,x2,t)  ( 0.075 * cos(0.000067597751162*t - pi/180*194.806 ) ...
                              + 0.095 * cos(0.000072921165921*t - pi/180*206.265 ) ...
                              + 0.10  * cos(0.000137879713787*t - pi/180*340.0   ) ...
                              + 0.395 * cos(0.000140518917083*t - pi/180*  0.0   ) ...
                              + 0.06  * cos(0.000145444119418*t - pi/180*42.97180) ) * (x1==x1);
    
    pd.configADCIRC.NBFR = 0;

    pd.isRivCont = true;
    n = [117.25; 503.5];
    n = n / norm(n);
    pd.xiRivCont = @(x1,x2,t) 0*x1;
    pd.uRivCont = @(x1,x2,t) -0.5 * n(1) * (x1==x1);
    pd.vRivCont = @(x1,x2,t) -0.5 * n(2) * (x1==x1);
  otherwise
    error('Invalid model.')
end % switch

end % function
