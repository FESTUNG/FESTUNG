function problemData = configureProblem(problemData)

%% Configuration to use: 
% - 'rotation' calls configureRotation()
% - 'analytical' calls configureAnalyticalTest()
problemData = setdefault(problemData, 'configSource', 'analytical');

% Polynomial approximation order
p = 0;

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
problemData.sweData.p = p;

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

problemData.transportData.p = p;
problemData.transportData.ordRK = problemData.sweData.schemeOrder; % as of now both models have to use the same RK method
problemData.transportData.tEnd = problemData.sweData.tEnd;
problemData.transportData.numSteps = problemData.sweData.numSteps; % as of now both models have to use the same time step
problemData.transportData.isVisGrid = false; % visualization of grid

problemData.transportData = execin('transport/configureProblem', problemData.transportData);
end % function