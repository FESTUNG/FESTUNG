function [ interface, name, g, zbAlg, fcAlg, gConst, NDTVAR, NOLIBF, NWP, bottomFric, F0, F0Alg, rhsAlg, F1Alg, F2Alg, xiOSAlg, t0, tEnd, dt, numSteps, ...
					 IRK, ISLOPE, ITRANS, CONVCR, minTol, output, isVisParam, isVisu, isVisuH, isVisv, isVisvH, isVisxi, isSolAvail, xi, u, v, tidalDomain, F, fT, NRAMP, ramping, ...
					 xiRI, uRI, vRI, rhsOSAlg, NBFR, xiOSX, xiOST, useStations, NSTAE, triE, coordE, NSTAV, triV, coordV, NHSTAR, NHSINC ] = userInput(problem, refinement)
switch problem
	case 0 % for debugging purposes
		interface = 'none';
    name      = 'debug';
    %% Domain specification
    % Must specify if an ADCIRC grid should be used
    gridGen   = 'manual';      
    coordSys  = 'unnecessary';
    % Must provide the underlying geometry and specify the bathymetry
    g = domainSquare(1);
		g.idE = zeros(g.numE,1);
		g.idE(g.baryE(:, 2) == 0) = 1; % south
		g.idE(g.baryE(:, 1) == 1) = 4; % east
		g.idE(g.baryE(:, 2) == 1) = 1; % north
		g.idE(g.baryE(:, 1) == 0) = 1; % west
		g.idE0T = g.idE(g.E0T);
%     zbAlg  = @(x1,x2) -0.1*(x1==x1);
		zbAlg  = @(x1,x2) -0.1*(1-x1<x2)-0.1;
    %% Coefficients, right hand sides, and boundary conditions
    fcAlg    = @(x1,x2  )  0*x1;
    gConst   = 9.81;
    bottomFric = 0;
    NOLIBF = 1; NWP = 0;
    
    % solution
    xi = @(x1,x2,t) 0*x1;
     u = @(x1,x2,t) 0*x1;
     v = @(x1,x2,t) 0*x1;
		xiOSAlg = @(x1,x2,t) (x1==x1);
    NRAMP = 0;
		xiRI = []; uRI = []; vRI = []; F1Alg = []; F2Alg = [];
    NDTVAR = 0; ITRANS = 0; CONVCR = []; NHSTAR = 0; NHSINC = [];
		IRK = 1; ISLOPE = [];
    
    F0    = 0;
    F0Alg = [];
    rhsAlg   = 0;
    %% Time stepping parameters
    t0       = 0;
    tEnd     = 0.01;
    numSteps = 5;
    %% Postprocessing parameter
    minTol      = 0.001;            % values of c1 in each vertex must be at least as high as this tolerance 
    %% Output parameters
    isVisGrid   = false;            % visualization of grid
		isVisParam  = false;						% visualization of zb and fc
    numPlots    = 20;               % number of plots for whole visualization
    isVisu      = false;            % visualization of x-velocity component
    isVisuH     = false;            % visualization of x-momentum component
    isVisv      = false;            % visualization of y-velocity component
    isVisvH     = false;            % visualization of y-momentum component
    isVisxi     = false;            % visualization of free surface elevation
    isSolAvail  = false;            % if solution available compute the error of the method
                                    % user must then create function handles H, u and v parametrizing the height and both velocity components
    %% station parameter
		useStations = false;						% whether or not to use elevation and velocity recording stations		
  case 1 % analytical test case
    %% General parameter
    % specifies usage of interface, 'none' means user has to specify all input manually
    interface = 'none';
    name      = 'analytical_test';
    %% Domain specification
    % Must specify if an ADCIRC grid should be used
    gridGen   = 'manual';      
    coordSys  = 'unnecessary';
    % Must provide the underlying geometry and specify the bathymetry
    X1 = [0 1 1 0]; X2 = [0 0 1 1];
		g = domainHierarchy(X1, X2, 0.3, refinement); % provides a nice unstructured mesh of maximal element size 0.3 of the unit square
    g.idE = zeros(g.numE,1);
		g.idE(g.baryE(:, 2) == 0) = 1; % south
		g.idE(g.baryE(:, 1) == 1) = 4; % east
		g.idE(g.baryE(:, 2) == 1) = 1; % north
		g.idE(g.baryE(:, 1) == 0) = 1; % west
    g.idE0T = g.idE(g.E0T);
    height = 0.05;
    A = 0.1;
    B = 0.1;
    C = 0.1;
    zbAlg  = @(x1,x2) -height*(2-x1-x2); % TODO non-zero 
    %% Coefficients, right hand sides, and boundary conditions
    fcAlg    = @(x1,x2  )  x1;
    gConst   = 9.81;
    bottomFric = 1;
    NOLIBF = 1; NWP = 0; NRAMP = 0;
    NDTVAR = 0; ITRANS = 0; CONVCR = []; NHSTAR = 0; NHSINC = [];
		IRK = 1; ISLOPE = [];
    
    % solution
    xi = @(x1,x2,t) C*(cos(0.5*pi*(x1-t)) + cos(0.5*pi*(x2-t))) - height*(2-x1-x2); % note that with this zb the value at point (1,1) is zero
     u = @(x1,x2,t) A*sin(pi*x1)*cos(2*pi*t);
     v = @(x1,x2,t) B*sin(pi*x2)*cos(2*pi*t);
    
    % auxilliary functions
    H   = @(x1,x2,t) xi(x1,x2,t) - zbAlg(x1,x2);
    u_t = @(x1,x2,t)   -2*A*pi*sin(pi*x1)*sin(2*pi*t);
    u_x = @(x1,x2,t)      A*pi*cos(pi*x1)*cos(2*pi*t);
    v_t = @(x1,x2,t)   -2*B*pi*sin(pi*x2)*sin(2*pi*t);
    v_y = @(x1,x2,t)      B*pi*cos(pi*x2)*cos(2*pi*t);
    H_t = @(x1,x2,t)  0.5*C*pi * ( sin(0.5*pi*(x1-t)) + sin(0.5*pi*(x2-t)) );
    H_x = @(x1,x2,t) -0.5*C*pi *   sin(0.5*pi*(x1-t))                       ;
    H_y = @(x1,x2,t) -0.5*C*pi *                       sin(0.5*pi*(x2-t))   ;
    
    F0    = 1;
    F0Alg = @(x1,x2,t) H_t(x1,x2,t) + (u_x(x1,x2,t) + v_y(x1,x2,t)) .* H(x1,x2,t) + u(x1,x2,t) .* H_x(x1,x2,t) + v(x1,x2,t) .* H_y(x1,x2,t);
    rhsAlg   = 1;
    F1Alg = @(x1,x2,t) u_t(x1,x2,t) .* H(x1,x2,t) + u(x1,x2,t) .* H_t(x1,x2,t) ...
                       + ( 2 * u(x1,x2,t) .* u_x(x1,x2,t) + gConst * H_x(x1,x2,t) ) .* H(x1,x2,t) ... 
                       + u(x1,x2,t) .* u(x1,x2,t) .* H_x(x1,x2,t) + u(x1,x2,t) .* v_y(x1,x2,t) .* H(x1,x2,t) ... 
                       + u(x1,x2,t) .* v(x1,x2,t) .* H_y(x1,x2,t) + gConst * height * H(x1,x2,t) ...
                       + bottomFric * ( u(x1,x2,t) .* u(x1,x2,t) + v(x1,x2,t) .* v(x1,x2,t) ).^0.5 .* u(x1,x2,t) ...
                       - fcAlg(x1,x2) .* v(x1,x2,t) .* H(x1,x2,t);
    F2Alg = @(x1,x2,t) v_t(x1,x2,t) .* H(x1,x2,t) + v(x1,x2,t) .* H_t(x1,x2,t) ...
                       + u_x(x1,x2,t) .* v(x1,x2,t) .* H(x1,x2,t) + u(x1,x2,t) .* v(x1,x2,t) .* H_x(x1,x2,t) + ...
                       + 2 * v(x1,x2,t) .* v_y(x1,x2,t) .* H(x1,x2,t) + v(x1,x2,t) .* v(x1,x2,t) .* H_y(x1,x2,t) ... 
                       + gConst * H_y(x1,x2,t) .* H(x1,x2,t) + gConst * height * H(x1,x2,t) ...
                       + bottomFric * ( u(x1,x2,t) .* u(x1,x2,t) + v(x1,x2,t) .* v(x1,x2,t) ).^0.5 .* v(x1,x2,t) ...
                       + fcAlg(x1,x2) .* u(x1,x2,t) .* H(x1,x2,t);
    %% Time stepping parameters
    t0       = 0;
    tEnd     = 1;
    numSteps = 150;
    %% Postprocessing parameter
    minTol      = 0.001;            % values of c1 in each vertex must be at least as high as this tolerance 
    %% Output parameters
    isVisGrid   = false;            % visualization of grid
		isVisParam  = false;						% visualization of zb and fc
    numPlots    = 20;               % number of plots for whole visualization
    isVisu      = true;             % visualization of x-velocity component
    isVisuH     = false;            % visualization of x-momentum component
    isVisv      = true;             % visualization of y-velocity component
    isVisvH     = false;            % visualization of y-momentum component
    isVisxi     = true;             % visualization of free surface elevation
    isSolAvail  = true;             % if solution available compute the error of the method
                                    % user must then create function handles H, u and v parametrizing the height and both velocity components
    %% station parameter
		useStations = false;						% whether or not to use elevation and velocity recording stations
  case 2 % manual input of bahamas parameters with only the grid constructed from ADCIRC files
    %% General parameter
    % specifies usage of interface, 'none' means user has to specify all input manually
    interface = 'none';           % Must provide type of grid generation, name, existence of rhs and if known type of used coordinate system
    name      = 'bahamas';
    %% Domain specification
    gridGen   = 'ADCIRCgrid';
    coordSys  = 'cartesian';
    %% Coefficients, right hand sides, and boundary conditions
    fcAlg    = @(x1,x2  ) 3.19 * 10^-5 * (x1==x1);
    gConst   = 9.81;
    bottomFric = 0.009;
    NOLIBF = 1; NWP = 0;
    F0       = 0;
    F0Alg    = @(x1,x2,t) zeros(size(x1));
    rhsAlg   = 0;
    F1Alg    = @(x1,x2,t) zeros(size(x1));
    F2Alg    = @(x1,x2,t) zeros(size(x1));
    NBFR     = 5;
    NRAMP    = 0;
    NDTVAR = 0; ITRANS = 0; CONVCR = []; NHSTAR = 0; NHSINC = [];
		IRK = 1; ISLOPE = [];
    xiOSAlg  = @(x1,x2,t) 0.075 * cos(0.000067597751162*t - pi/180*194.806 ) ...
                        + 0.095 * cos(0.000072921165921*t - pi/180*206.265 ) ...
                        + 0.10  * cos(0.000137879713787*t - pi/180*340.0   ) ...
                        + 0.395 * cos(0.000140518917083*t - pi/180*  0.0   ) ...
                        + 0.06  * cos(0.000145444119418*t - pi/180*42.97180);
    %% Time stepping parameters
    t0       = 0;
    tEnd     = 12*24*3600;
    dt       = 15;
    %% Postprocessing parameter
    minTol      = 0.001;            % values of c1 in each vertex must be at least as high as this tolerance 
    %% Output parameters
    isVisGrid   = false;            % visualization of grid
		isVisParam  = false;						% visualization of zb and fc
    numPlots    = 240;              % number of plots for whole visualization
    isVisu      = true;             % visualization of x-velocity component
    isVisuH     = false;            % visualization of x-momentum component
    isVisv      = true;             % visualization of y-velocity component
    isVisvH     = false;            % visualization of y-momentum component
    isVisxi     = true;             % visualization of free surface elevation
    isSolAvail  = false;            % if solution available compute the error of the method
                                    % user must then create function handles H, u and v parametrizing the height and both velocity components
		%% station parameter
		useStations = true; 						% whether or not to use elevation and velocity recording stations
		NSTAE = 4;          						% here we specify stations manually
		XEL = [38666.66 56097.79 41262.60 59594.66];
    YEL = [49333.32  9612.94 29775.73 41149.62];
    NSTAV = NSTAE;
		XEV = XEL;
    YEV = YEL;
  case 3 % arbitrary ADCIRC input
    %% General parameters
    % specifies usage of interface 'ADCIRC' means everything is provided from ADCIRC parameter files only name of files has to be specified
    interface = 'ADCIRC';
    name      = 'gom3k';
    %% Domain specification
    gridGen   = 'unnecessary';
    coordSys  = 'unnecessary';
    %% Coefficients, right hand sides, and boundary conditions
    F1Alg   = [];
    F2Alg   = [];
    xiOSAlg = [];
    %% Output parameters
    isVisGrid   = false;            % visualization of grid
		isVisParam  = false;						% visualization of zb and fc
    numPlots    = 200;              % number of plots for whole visualization
    isVisu      = true;             % visualization of x-velocity component
    isVisuH     = false;            % visualization of x-momentum component
    isVisv      = true;             % visualization of y-velocity component
    isVisvH     = false;            % visualization of y-momentum component
    isVisxi     = true;             % visualization of free surface elevation
    isSolAvail  = false;            % if solution available compute the error of the method
                                    % user must then create function handles H, u and v parametrizing the height and both velocity components
		%% station parameter
		useStations = false;						% whether or not to use elevation and velocity recording stations
  case 4 % arbitrary ADCIRC input
    %% General parameters
    % specifies usage of interface 'ADCIRC' means everything is provided from ADCIRC parameter files only name of files has to be specified
    interface = 'ADCIRC';
    name      = 'test2';
    %% Domain specification
    gridGen   = 'unnecessary';
    coordSys  = 'unnecessary';
    %% Coefficients, right hand sides, and boundary conditions
    F1Alg   = [];
    F2Alg   = [];
    xiOSAlg = [];
    %% Output parameters
    isVisGrid   = false;            % visualization of grid
		isVisParam  = false;						% visualization of zb and fc
    numPlots    = 72;               % number of plots for whole visualization
    isVisu      = true;             % visualization of x-velocity component
    isVisuH     = false;            % visualization of x-momentum component
    isVisv      = true;             % visualization of y-velocity component
    isVisvH     = false;            % visualization of y-momentum component
    isVisxi     = true;             % visualization of free surface elevation
    isSolAvail  = false;            % if solution available compute the error of the method
                                    % user must then create function handles H, u and v parametrizing the height and both velocity components
		%% station parameter
		useStations = false;						% whether or not to use elevation and velocity recording stations
  case 5 % arbitrary ADCIRC input
    %% General parameters
    % specifies usage of interface 'ADCIRC' means everything is provided from ADCIRC parameter files only name of files has to be specified
    interface = 'ADCIRC';
    name      = 'gom';
    %% Domain specification
    gridGen   = 'unnecessary';
    coordSys  = 'unnecessary';
    %% Coefficients, right hand sides, and boundary conditions
    F1Alg   = [];
    F2Alg   = [];
    xiOSAlg = [];
    %% Output parameters
    isVisGrid   = false;            % visualization of grid
		isVisParam  = false;						% visualization of zb and fc
    numPlots    = 200;              % number of plots for whole visualization
    isVisu      = true;             % visualization of x-velocity component
    isVisuH     = false;            % visualization of x-momentum component
    isVisv      = true;             % visualization of y-velocity component
    isVisvH     = false;            % visualization of y-momentum component
    isVisxi     = true;             % visualization of free surface elevation
    isSolAvail  = false;            % if solution available compute the error of the method
                                    % user must then create function handles H, u and v parametrizing the height and both velocity components
		%% station parameter
		useStations = false;						% whether or not to use elevation and velocity recording stations
  case 6
    %% General parameters
    % specifies usage of interface 'ADCIRC' means everything is provided from ADCIRC parameter files only name of files has to be specified
    interface = 'ADCIRC';
    name      = 'east';
    %% Domain specification
    gridGen   = 'unnecessary';
    coordSys  = 'unnecessary';
    %% Coefficients, right hand sides, and boundary conditions
    F1Alg   = [];
    F2Alg   = [];
    xiOSAlg = [];
    %% Output parameters
    isVisGrid   = false;            % visualization of grid
		isVisParam  = false;						% visualization of zb and fc
    numPlots    = 200;              % number of plots for whole visualization
    isVisu      = true;             % visualization of x-velocity component
    isVisuH     = false;            % visualization of x-momentum component
    isVisv      = true;             % visualization of y-velocity component
    isVisvH     = false;            % visualization of y-momentum component
    isVisxi     = true;             % visualization of free surface elevation
    isSolAvail  = false;            % if solution available compute the error of the method
                                    % user must then create function handles H, u and v parametrizing the height and both velocity components
		%% station parameter
		useStations = false;						% whether or not to use elevation and velocity recording stations
  otherwise
    error('Unknown problem.');
end % switch

%% Do not change below!
fprintf('\n SOLVING PROBLEM: %s \n\n', upper(name));
%% create folders for output
warning('off', 'all');
mkdir(['output_' name]);
warning('on', 'all');

%% grid construction
fprintf('Construct grid and read parameters.\n');
if strcmp(interface, 'ADCIRC')
  [ g, fcAlg, zbAlg, NOLIBF, NWP, bottomFric, tidalDomain, NRAMP, F, fT, ramping, xiRI, uRI, vRI, NBFR, xiOSX, xiOST, gConst, NDTVAR, dt, t0, tEnd, ...
		IRK, ISLOPE, ITRANS, CONVCR, minTol, NSTAE, XEL, YEL, NSTAV, XEV, YEV, NHSTAR, NHSINC ] = fort2Mat(name);
  rhsOSAlg = 0;
elseif strcmp(interface, 'none')
  tidalDomain = 0;
  rhsOSAlg = 1;
  F = []; fT = []; xiOSX = []; xiOST = []; 
  if strcmp(gridGen, 'manual')
    NBFR = [];
  elseif strcmp(gridGen, 'ADCIRCgrid')
    if strcmp(coordSys, 'unknown')
      [ICS, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, SLAM0, SFEA0] = fort15Read(['fort_' name '.15']);
      [g, DP] = fort14Read(['fort_' name '.14'], ICS, SLAM0, SFEA0);
    elseif strcmp(coordSys, 'cartesian')
      [g, DP] = fort14Read(['fort_' name '.14'], 1);
    else
      error('Invalid coordinate sytem specification.');
    end % if
    [g, NEDGES, NEDNO, NLEDN, NRAEDN, NRIEDN, ETRI, UNRI, UTRI, NSEDN] = fort17Read(['fort_' name '.17'], g, NBFR);
    [g, NLEDN, NRAEDN, NRIEDN, NSEDN] = generateGridDataFromADCIRC(g, NEDGES, NEDNO, NLEDN, NRAEDN, NRIEDN, NSEDN);
		g.idE = zeros(g.numE, 1);
		g.idE(NLEDN ) = 1;
		g.idE(NRAEDN) = 2;
		g.idE(NRIEDN) = 3;
		g.idE(NSEDN ) = 4;
		g.idE0T = g.idE(g.E0T); % local edge IDs
		xiRI = zeros(g.numT, 1);
		 uRI = zeros(g.numT, 1);
		 vRI = zeros(g.numT, 1);
		xiRI(g.T0E(NRIEDN,1)) = ETRI;
		 uRI(g.T0E(NRIEDN,1)) = UNRI .* g.nuE(NRIEDN, 1) - UTRI .* g.nuE(NRIEDN, 2);
		 vRI(g.T0E(NRIEDN,1)) = UNRI .* g.nuE(NRIEDN, 2) + UTRI .* g.nuE(NRIEDN, 1);
		if ~isequal(g.T0E(NSEDN,2), zeros(length(NSEDN),1))
      error('The marked open sea edges are not located on the boundary.');
    end % if
    zbAlg = @(x1,x2) evaluateFuncFromVertexValues(g,-DP,x1,x2);
		assert( max( max( abs( DP(g.V0T) + zbAlg(g.coordV0T(:,:,1), g.coordV0T(:,:,2)) ) ) ) < 10^-5, 'Bathymetry not constructed correctly.' ); % unusal large tol
  else
    error('Invalid grid construction option.');
  end % if
else 
  error('Invalid interface.');
end % if
if isVisGrid
  visualizeGrid(g);
end % if

%% time stepping specification
fprintf('Determine length and number of timesteps for simulation.\n')
if     exist('dt', 'var') && ~exist('numSteps', 'var')
  numSteps = round((tEnd-t0) / dt);
elseif ~exist('dt', 'var') && exist('numSteps', 'var')
  dt = (tEnd-t0) / numSteps;
elseif  exist('dt', 'var') && exist('numSteps', 'var')
  if ((tEnd-t0) / numSteps ~= dt)
    dt = (tEnd-t0) / numSteps;
    warning('dt has been reset as %g to account for length of simulation and number of time steps.', dt);
  end % if
else
  error('Either time increment or number of time steps has to be specified.');
end % if
numPlots = min(numPlots, numSteps);
output = max(floor(round(numSteps / numPlots * 1000) / 1000), 1); % for plotting only after certain times
if output ~= numSteps / numPlots
	warning('The number of created plots could differ from the specified number of plots.');
end % if

if isSolAvail
  if exist('xiOSAlg', 'var')
    warning('Specification of open sea boundary condition for free surface elevation might be in conflict with analytical solution. Please check!');
  else 
    % boundary condition specification
    xiOSAlg = @(x1,x2,t) xi(x1,x2,t);
  end % if
	if exist('xiRIAlg', 'var')
    warning('Specification of river boundary condition for free surface elevation might be in conflict with analytical solution. Please check!');
  else 
    % boundary condition specification
    xiRI = @(x1,x2,t) xi(x1,x2,t);
  end % if
	if exist('uRIAlg', 'var')
    warning('Specification of river boundary condition for x-velocity component might be in conflict with analytical solution. Please check!');
  else 
    % boundary condition specification
     uRI = @(x1,x2,t)  u(x1,x2,t);
  end % if
	if exist('vRIAlg', 'var')
		warning('Specification of river boundary condition for y-velocity component might be in conflict with analytical solution. Please check!');
  else 
    % boundary condition specification
     vRI = @(x1,x2,t)  v(x1,x2,t); % TODO specification of river boundary with function handles
  end % if
else
  xi = []; u = []; v = [];
  % no right hand side for continuity equation
  F0 = 0;
  F0Alg = [];
end % if

if ~exist('rhsAlg', 'var')
  rhsAlg = 0;
elseif strcmp(interface, 'ADCIRC')
  if rhsAlg == 1
    warning('Possible conflicts for right hand sides of momentum equations. Algebraic right hand sides are being neglected.');
  end % if
  rhsAlg = 0;
end % if

% various
if ~NRAMP
  ramping = [];
end % if

%% stations
if useStations
	if NDTVAR == 1
		error('For variable time step size the use of output stations is not supported.');
	end % if
	if NSTAE == 0 && NSTAV == 0
		warning('No stations are specified.')
		useStations = false;
	else
		coordE = [XEL; YEL].';
		triE = cell(1, NSTAE);
		for i = 1:NSTAE
			triE{i} = findTrinagle(g, coordE(i,1), coordE(i,2));
		end % for
		coordV = [XEV; YEV].';
		triE = cell2mat(triE);
		triV = cell(1, NSTAV);
		for i = 1:NSTAV
			triV{i} = findTrinagle(g, coordV(i,1), coordV(i,2));
		end % for
		triV = cell2mat(triV);
	end % if
end % if
if ~useStations
  NSTAE = [];
  NSTAV = [];
  triE = [];
  coordE = [];
  triV = [];
  coordV = [];
end % if
end % function
