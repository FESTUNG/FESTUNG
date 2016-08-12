function problemData = configureProblem(problemData)

%% Configuration to use: 
% - 'rotation' calls configureRotation()
% - 'analytical' calls configureAnalyticalTest()
problemData = setdefault(problemData, 'configSource', 'analytical');

%% What kind of grid to use:
% - 'square' creates a unit square [0,1]x[0,1] with given pd.hmax,
%   open sea boundary in the east (type 4), and land boundary (type 1) on 
%   all other edges 
% - 'hierarchical' creates a unit square [0,1]x[0,1] with specified hmax
%   and performs uniform refinement according to parameter 'refinement'.
%   Boundary type 4 on east-boundary, 1 on all others.
% - 'ADCIRC' reads grid information from 'swe/fort_<name>.{14,17}'.
problemData = setdefault(problemData, 'gridSource', 'hierarchical');
problemData = setdefault(problemData, 'refinement', 3);
% problemData = setdefault(problemData, 'hmax', 0.3);

% Polynomial approximation order
problemData = setdefault(problemData, 'p', 1);

% Runge-Kutta order
problemData = setdefault(problemData, 'ordRK', min(problemData.p+1,3)); % as of now both models have to use the same RK method

% Time stepping specification, as of now both models have to use the same number of time steps
switch problemData.configSource
  case 'rotation'
    problemData = setdefault(problemData, 'hmax', 2^-6);
    problemData = setdefault(problemData, 'tEnd', (100/3142)*2*pi);
    problemData = setdefault(problemData, 'numSteps', 100);
  case 'analytical'
    problemData = setdefault(problemData, 'hmax', 1);
    problemData = setdefault(problemData, 'tEnd', 5000);
    problemData = setdefault(problemData, 'numSteps', 10*2^problemData.refinement*(problemData.p+1));
  otherwise
    error('Invalid config source.')
end % switch

% Configuration for shallow water solver
problemData.sweData = struct;

% Specification of configuration type for shallow water model
switch problemData.configSource
  case 'rotation'
    problemData.sweData.configSource = 'debug';
  case 'analytical'
    problemData.sweData.configSource = 'analytical';
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

problemData.sweData = execin('swe/configureProblem', problemData.sweData);

% Configuration for transport solver
problemData.transportData = struct;

% Specification of configuration type for transport model and solution
switch problemData.configSource
  case 'rotation'
    problemData.transportData.configSource = 'rotation';
    
  case 'analytical'
    problemData.transportData.configSource = 'analytical';
    
    % analytical functions from swe necessary for convergence test
    problemData.transportData.hCont = problemData.sweData.hCont;
    problemData.transportData.h_tCont = problemData.sweData.h_tCont;
    problemData.transportData.h_xCont = problemData.sweData.h_xCont;
    problemData.transportData.h_yCont = problemData.sweData.h_yCont;
    problemData.transportData.uCont = problemData.sweData.uCont;
    problemData.transportData.u_xCont = problemData.sweData.u_xCont;
    problemData.transportData.vCont = problemData.sweData.vCont;
    problemData.transportData.v_yCont = problemData.sweData.v_yCont;
  otherwise
    error('Invalid config source.')
end % switch

problemData.transportData.gridSource = problemData.gridSource;
problemData.transportData.refinement = problemData.refinement;
problemData.transportData.hmax = problemData.hmax;
problemData.transportData.p = problemData.p;

problemData.transportData.ordRK = problemData.ordRK;
problemData.transportData.tEnd = problemData.tEnd;
problemData.transportData.numSteps = problemData.numSteps;

problemData.transportData.isVisGrid = false; % visualization of grid

problemData.transportData = execin('transport/configureProblem', problemData.transportData);
end % function