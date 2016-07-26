
function pd = configureProblem(pd)
%% Name of the problem
pd.name = 'inf'; 

%% Polynomial approximation order
% Piecewise constant (0), piecewise linear (1), or piecewise quadratic (2)
% p finite direction
% pInf infinite direction
% pInf >= p
pd.p = 2;
pd.pInf = 2;


%% Grid Parameter, beta = 1/L where L is a typical length scale of the problem
pd.beta = 0.01;

%% Time stepping parameters
pd.schemeType = 'explicit'; % type of time stepping scheme ('explicit')

%% Model parameters
% Some may be overwritten by fort.15 config files
pd.typeFlux = 'Lax-Friedrichs'; % Type of interior flux ('Lax-Friedrichs', 'Roe')
pd.isRiemOS = true; % Riemann solver type on open sea boundary ('Lax-Friedrichs', 'Roe', or 'none')
pd.typeBdrL = 'riemann'; % Flux type on land boundary ('reflected', 'natural', or 'riemann')

%% Testcase
pd = configureAnalyticalTest2(pd);

%% Waitbar
pd.isWaitbar = false;

end % function



%% Analytical solution
function pd = configureAnalyticalTest(pd)
pd.isSolutionAvail = true;
pd.isRhsAvail = true;
pd.isTidalDomain = false;

% Overwrite grid parameters
pd.gridSource = 'hierarchical';
pd.isSpherical = false; 
pd.hmax = 30; % Maximum element size of initial grid
pd.refinement = 0;  % Grid refinement level

% Overwrite time-stepping parameters
pd.t0 = 0; % Start time of simulation
pd.tEnd = 1; % End time of simulation
pd.numSteps = 1; % Number of time steps

pd.isAdaptiveTimestep = false; % Use adaptive timestep width
pd.dt = (pd.tEnd - pd.t0) / pd.numSteps;

pd.isSteadyState = false; % End simulation upon convergence

% Solution parameters
height = 1;
A = 0.01;
pd.gConst = 9.81;

pd.isBottomFrictionNonlinear = true; % NOLIBF
pd.isBottomFrictionVarying = false; % NWP
pd.bottomFrictionCoef = 1;

% Ramping function, bathymetry, and Coriolis coefficient
pd.isRamp = false;
pd.zbCont = @(x1,x2) -height * 10 *ones(size(x2));
pd.fcCont = @(x1,x2) zeros(size(x1));

% Analytical solution
pd.xiCont = @(x1,x2,t) zeros(size(x1));
pd.uCont = @(x1,x2,t) A*x1;
pd.vCont = @(x1,x2,t) A*A*x2.*(x2-100);

% Auxiliary functions (derivatives etc.)
pd.hCont = @(x1,x2,t) pd.xiCont(x1,x2,t) - pd.zbCont(x1,x2);
pd.u_tCont  = @(x1,x2,t)    zeros(size(x1));
pd.u_xCont  = @(x1,x2,t)   A*ones(size(x1));
pd.u_yCont  = @(x1,x2,t)    zeros(size(x1));
pd.v_tCont  = @(x1,x2,t)    zeros(size(x1));
pd.v_xCont  = @(x1,x2,t)    zeros(size(x1));
pd.v_yCont  = @(x1,x2,t)    A*A*(x2-100+x2);
pd.zb_xCont = @(x1,x2,t)    zeros(size(x1));
pd.zb_yCont = @(x1,x2,t)    zeros(size(x1));
pd.h_tCont  = @(x1,x2,t)    zeros(size(x1));
pd.h_xCont  = @(x1,x2,t)    zeros(size(x1));
pd.h_yCont  = @(x1,x2,t)    zeros(size(x1));

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
                       pd.u_yCont(x1,x2,t) .* pd.vCont(x1,x2,t) .* pd.hCont(x1,x2,t) + ...
                       pd.gConst * pd.zb_xCont(x1,x2,t) .* pd.hCont(x1,x2,t) + ...
                       pd.bottomFrictionCoef * sqrt( pd.uCont(x1,x2,t) .* pd.uCont(x1,x2,t) + pd.vCont(x1,x2,t) .* pd.vCont(x1,x2,t) ) .* pd.uCont(x1,x2,t) - ...
                       pd.fcCont(x1,x2) .* pd.vCont(x1,x2,t) .* pd.hCont(x1,x2,t);
pd.f2Cont = @(x1,x2,t) pd.v_tCont(x1,x2,t) .* pd.hCont(x1,x2,t) + ...
                       pd.vCont(x1,x2,t) .* pd.h_tCont(x1,x2,t) + ...
                       pd.u_xCont(x1,x2,t) .* pd.vCont(x1,x2,t) .* pd.hCont(x1,x2,t) + ...
                       pd.uCont(x1,x2,t) .* pd.vCont(x1,x2,t) .* pd.h_xCont(x1,x2,t) + ...
                       pd.uCont(x1,x2,t) .* pd.v_xCont(x1,x2,t) .* pd.hCont(x1,x2,t) + ...
                       2 * pd.vCont(x1,x2,t) .* pd.v_yCont(x1,x2,t) .* pd.hCont(x1,x2,t) + ...
                       pd.vCont(x1,x2,t) .* pd.vCont(x1,x2,t) .* pd.h_yCont(x1,x2,t) + ... 
                       pd.gConst * pd.h_yCont(x1,x2,t) .* pd.hCont(x1,x2,t) + ...
                       pd.gConst * pd.zb_yCont(x1,x2,t) .* pd.hCont(x1,x2,t) + ...
                       pd.bottomFrictionCoef * sqrt( pd.uCont(x1,x2,t) .* pd.uCont(x1,x2,t) + pd.vCont(x1,x2,t) .* pd.vCont(x1,x2,t) ) .* pd.vCont(x1,x2,t) + ...
                       pd.fcCont(x1,x2) .* pd.uCont(x1,x2,t) .* pd.hCont(x1,x2,t);
                    
                     
% Boundary conditions
pd.xiOSCont = pd.xiCont;
pd.isRivCont = false;

%zb is linear in finite and constant in infinite direction
%OS-Bdr at x1 = 100
pd.zbOSAlg = @(x1,x2) pd.zbCont(100,x2);
pd.HOSAlg  = @(x1,x2,t) pd.xiCont(x1,x2,t) - pd.zbOSAlg(x1,x2);

%height >= 0
pd.minTol = 0.0001;
end % function

function pd = configureAnalyticalTest2(pd)
pd.isSolutionAvail = true;
pd.isRhsAvail = true;
pd.isTidalDomain = false;

% Overwrite grid parameters
pd.gridSource = 'hierarchical';
pd.isSpherical = false; 
pd.hmax = 30; % Maximum element size of initial grid
pd.refinement = 0;  % Grid refinement level

% Overwrite time-stepping parameters
pd.t0 = 0; % Start time of simulation
pd.tEnd = 1; % End time of simulation
pd.numSteps = 1; % Number of time steps

pd.isAdaptiveTimestep = false; % Use adaptive timestep width
pd.dt = (pd.tEnd - pd.t0) / pd.numSteps;

pd.isSteadyState = false; % End simulation upon convergence

% Solution parameters
height = 1;
A = 0.01;
B = 0.1;
pd.gConst = 9.81;

pd.isBottomFrictionNonlinear = true; % NOLIBF
pd.isBottomFrictionVarying = false; % NWP
pd.bottomFrictionCoef = 1;

% Ramping function, bathymetry, and Coriolis coefficient
pd.isRamp = false;
pd.zbCont = @(x1,x2) -height * (6 - A * x2);
pd.fcCont = @(x1,x2) zeros(size(x1));

% Analytical solution
pd.xiCont = @(x1,x2,t) B * ( cos(A*pi*(x2)) );
pd.uCont = @(x1,x2,t) A*x1;
pd.vCont = @(x1,x2,t) B*sin(A*pi*x2);

% Auxiliary functions (derivatives etc.)
pd.hCont = @(x1,x2,t) pd.xiCont(x1,x2,t) - pd.zbCont(x1,x2);
pd.u_tCont  = @(x1,x2,t)    zeros(size(x1));
pd.u_xCont  = @(x1,x2,t)   A*ones(size(x1));
pd.u_yCont  = @(x1,x2,t)    zeros(size(x1));
pd.v_tCont  = @(x1,x2,t)    zeros(size(x1));
pd.v_xCont  = @(x1,x2,t)    zeros(size(x1));
pd.v_yCont  = @(x1,x2,t)    A*B*pi*cos(A*pi*x2);
pd.zb_xCont = @(x1,x2,t)    zeros(size(x1));
pd.zb_yCont = @(x1,x2,t)    A * height*ones(size(x2));
pd.h_tCont  = @(x1,x2,t)    zeros(size(x1));
pd.h_xCont  = @(x1,x2,t)    zeros(size(x1));
pd.h_yCont  = @(x1,x2,t)    -A*B*pi * sin(A*pi*(x2)) - A*height*ones(size(x2));

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
                       pd.u_yCont(x1,x2,t) .* pd.vCont(x1,x2,t) .* pd.hCont(x1,x2,t) + ...
                       pd.gConst * pd.zb_xCont(x1,x2,t) .* pd.hCont(x1,x2,t) + ...
                       pd.bottomFrictionCoef * sqrt( pd.uCont(x1,x2,t) .* pd.uCont(x1,x2,t) + pd.vCont(x1,x2,t) .* pd.vCont(x1,x2,t) ) .* pd.uCont(x1,x2,t) - ...
                       pd.fcCont(x1,x2) .* pd.vCont(x1,x2,t) .* pd.hCont(x1,x2,t);
pd.f2Cont = @(x1,x2,t) pd.v_tCont(x1,x2,t) .* pd.hCont(x1,x2,t) + ...
                       pd.vCont(x1,x2,t) .* pd.h_tCont(x1,x2,t) + ...
                       pd.u_xCont(x1,x2,t) .* pd.vCont(x1,x2,t) .* pd.hCont(x1,x2,t) + ...
                       pd.uCont(x1,x2,t) .* pd.vCont(x1,x2,t) .* pd.h_xCont(x1,x2,t) + ...
                       pd.uCont(x1,x2,t) .* pd.v_xCont(x1,x2,t) .* pd.hCont(x1,x2,t) + ...
                       2 * pd.vCont(x1,x2,t) .* pd.v_yCont(x1,x2,t) .* pd.hCont(x1,x2,t) + ...
                       pd.vCont(x1,x2,t) .* pd.vCont(x1,x2,t) .* pd.h_yCont(x1,x2,t) + ... 
                       pd.gConst * pd.h_yCont(x1,x2,t) .* pd.hCont(x1,x2,t) + ...
                       pd.gConst * pd.zb_yCont(x1,x2,t) .* pd.hCont(x1,x2,t) + ...
                       pd.bottomFrictionCoef * sqrt( pd.uCont(x1,x2,t) .* pd.uCont(x1,x2,t) + pd.vCont(x1,x2,t) .* pd.vCont(x1,x2,t) ) .* pd.vCont(x1,x2,t) + ...
                       pd.fcCont(x1,x2) .* pd.uCont(x1,x2,t) .* pd.hCont(x1,x2,t);
                    
                     
% Boundary conditions
pd.xiOSCont = pd.xiCont;
pd.isRivCont = false;

%zb is linear in finite and constant in infinite direction
%OS-Bdr at x1 = 100
pd.zbOSAlg = @(x1,x2) pd.zbCont(100,x2);
pd.HOSAlg  = @(x1,x2,t) pd.xiCont(x1,x2,t) - pd.zbOSAlg(x1,x2);

%height >= 0
pd.minTol = 0.0001;
end % function
