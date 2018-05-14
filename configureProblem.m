function problemData = configureProblem(problemData)
%% Parameters.
% Name of testcase
problemData = setdefault(problemData, 'testcase', 'debug_jump');

problemData.eta = 1; % penalty parameter (eta>0)

% Number of elements in x- and y-direction
problemData = setdefault(problemData, 'numElem', [50, 10]);

% Local polynomial approximation order (0 to 5)
problemData = setdefault(problemData, 'p', 0);

% Order of quadrature rule
problemData = setdefault(problemData, 'qOrd', 2*problemData.p + 1);
problemData = setdefault(problemData, 'qOrdMax', problemData.qOrd);

% Time stepping parameters
problemData = setdefault(problemData, 't0', 0);  % start time
problemData = setdefault(problemData, 'tEnd', 10000);  % end time
problemData = setdefault(problemData, 'numSteps', ceil(problemData.tEnd/0.5));  % number of time steps

% Discard time derivative and compute stationary solution
problemData = setdefault(problemData, 'isStationary', false);  
problemData = setdefault(problemData, 'isCoupling', false);   
problemData = setdefault(problemData, 'isJumpCoupling', true);

% Visualization settings
problemData = setdefault(problemData, 'isVisGrid', false);  % visualization of grid
problemData = setdefault(problemData, 'isVisSol', true);  % visualization of solution
problemData = setdefault(problemData, 'outputFrequency', 100); % no visualization of every timestep
problemData = setdefault(problemData, 'outputBasename', ...  % Basename of output files
                         ['output' filesep problemData.problemName '_' problemData.testcase ]); 
problemData = setdefault(problemData, 'outputTypes', { 'vtk' });  % Type of visualization files ('vtk, 'tec')

%% Parameter check.
assert(problemData.p >= 0 && problemData.p <= 5, 'Polynomial order must be zero to four.')
assert(problemData.numSteps > 0, 'Number of time steps must be positive.')

%% Coefficients and boundary data.
isHotstart = false;

switch problemData.testcase
  case 'debug_jump'
    isHotstart = false;
    hotstartFile = ['darcy_2dv' filesep 'debug_jump_p0_50x10.mat'];
    
    % width and height of computational domain
    zPMCont = @(x) -6 * ones(size(x));
    zBotCont = @(x) 0.75 * (cos((x-50)/30 * pi) + 1) .* (20 <= x & x <= 80);

    domainWidth = linspace(0, 100, problemData.numElem(1)+1);
    domainHeight = [zPMCont(domainWidth); zBotCont(domainWidth)];
    idDirichlet = 3; idNeumann = [1, 2, 4];
    
    k = 0.001;
    
    problemData.hDCont = @(t,x,z) 5 * ones(size(x));
    problemData.gNCont = @(t,x,z) zeros(size(x));
    problemData.fCont = @(t,x,z) zeros(size(x));
    problemData.KCont = cellfun(@(c) @(t,x,z) c * ones(size(x)), {k, 0; 0, k}, 'UniformOutput', false);
    
  case 'showcase'
    isHotstart = false;
    hotstartFile = ['darcy_2dv' filesep 'showcase_p0_50x10.mat'];
    
    % width and height of computational domain
    zPMCont = @(x) -4 * ones(size(x));
    zBotCont = @(x) (cos((x-50)/30 * pi) + 1) .* (20 <= x & x <= 80);

    domainWidth = linspace(0, 100, problemData.numElem(1)+1);
    domainHeight = [zPMCont(domainWidth); zBotCont(domainWidth)];
    idDirichlet = 3; idNeumann = [1, 2, 4];
    
    k = 0.001;
    
    problemData.h0Cont = @(x,z) 5 * ones(size(x));
    problemData.q10Cont = @(x,z) zeros(size(x));
    problemData.q20Cont = @(x,z) zeros(size(x));
    
    problemData.hDCont = @(t,x,z) 5 * ones(size(x));
    problemData.gNCont = @(t,x,z) zeros(size(x));
    problemData.fCont = @(t,x,z) zeros(size(x));
    problemData.KCont = cellfun(@(c) @(t,x,z) c * ones(size(x)), {k, 0; 0, k}, 'UniformOutput', false);
    
  case {'coupled_constXi', 'coupled_stationary', 'coupled_transient'}
    % Parameters
    aConst = 0;
    bConst = 0.005;
    betaConst = 0.3;
    etaConst = 0.003;
    rhoConst = 0.08;
    tauConst = 0.08;
    kConst = 0.01;
    kappaConst = 1;
    lambdaConst = 0.07;
    muConst = 0.07;
    
    % Width and height of computational domain
    zPMCont = @(x) -5 * ones(size(x));
    zBotCont = @(x) aConst + bConst * x;
    
    domainWidth = linspace(0, 100, problemData.numElem(1)+1);
    domainHeight = [zPMCont(domainWidth); zBotCont(domainWidth)];
    idDirichlet = [1, 2, 3, 4]; idNeumann = -1;
    
    % Analytical solution
    switch problemData.testcase
      case 'coupled_constXi'
        xiCont = @(t,x) 5 * ones(size(x));

        dxXiCont = @(t,x) zeros(size(x));
        dxdxXiCont = @(t,x) zeros(size(x));
        dtXiCont = @(t,x) zeros(size(x));
      case 'coupled_stationary'
        xiCont = @(t,x) 5 + etaConst * sin(rhoConst * x);

        dxXiCont = @(t,x) etaConst * rhoConst * cos(rhoConst * x);
        dxdxXiCont = @(t,x) -etaConst * rhoConst * rhoConst * sin(rhoConst * x);
        dtXiCont = @(t,x) zeros(size(x));
      case 'coupled_transient'
        xiCont = @(t,x) 5 + etaConst * sin(rhoConst * x + tauConst * t);

        dxXiCont = @(t,x) etaConst * rhoConst * cos(rhoConst * x + tauConst * t);
        dxdxXiCont = @(t,x) -etaConst * rhoConst * rhoConst * sin(rhoConst * x + tauConst * t);
        dtXiCont = @(t,x) etaConst * tauConst * cos(rhoConst * x + tauConst * t);
    end % switch
    
    switch problemData.testcase
      case {'coupled_constXi', 'coupled_stationary'}
        omegaCont = @(t,x) kappaConst * cos(lambdaConst * x);
        dtOmegaCont = @(t,x) zeros(size(x));
        dxOmegaCont = @(t,x) -lambdaConst * kappaConst * sin(lambdaConst * x);
        dxdxOmegaCont = @(t,x) - lambdaConst * lambdaConst * kappaConst * cos(lambdaConst * x);
      case {'coupled_transient'}
        omegaCont = @(t,x) kappaConst * cos(lambdaConst * x + muConst * t);
        dtOmegaCont = @(t,x) -muConst * kappaConst * sin(lambdaConst * x + muConst * t);
        dxOmegaCont = @(t,x) -lambdaConst * kappaConst * sin(lambdaConst * x + muConst * t);
        dxdxOmegaCont = @(t,x) -lambdaConst * lambdaConst * kappaConst * cos(lambdaConst * x + muConst * t);
    end % switch
    
    problemData.hCont = @(t,x,z) xiCont(t,x) + (sin(betaConst * z) - sin(betaConst * zBotCont(x))) .* omegaCont(t,x);
    problemData.q1Cont = @(t,x,z) -dxXiCont(t,x) ...
      + betaConst * bConst * cos(betaConst * zBotCont(x)) .* omegaCont(t,x) ...
      - (sin(betaConst * z) - sin(betaConst * zBotCont(x))) .* dxOmegaCont(t,x);
    problemData.q2Cont = @(t,x,z) -betaConst * cos(betaConst * z) .* omegaCont(t, x);
    
    % Diffusion matrix
    problemData.KCont = cellfun(@(c) @(t,x,z) c * ones(size(x)), {kConst, 0; 0, kConst}, 'UniformOutput', false);
    
    % Derivatives
    dThCont = @(t,x,z) dtXiCont(t,x) + (sin(betaConst * z) - sin(betaConst * zBotCont(x))) .* dtOmegaCont(t,x);
    dXhCont = @(t,x,z) -problemData.q1Cont(t,x,z);
    dZhCont = @(t,x,z) -problemData.q2Cont(t,x,z);
    dXdXhCont = @(t,x,z) -dxdxXiCont(t,x) ...
      + betaConst * betaConst * bConst * bConst * sin(betaConst * zBotCont(x)) .* omegaCont(t,x) ...
      - 2 * betaConst * bConst * cos(betaConst * zBotCont(x)) .* dxOmegaCont(t,x) ...
      + (sin(betaConst * z) - sin(betaConst * zBotCont(x))) .* dxdxOmegaCont(t,x);
    dZdZhCont = @(t,x,z) -betaConst * betaConst * sin(betaConst * z) .* omegaCont(t, x);
    dXdZhCont = @(t,x,z) betaConst * cos(betaConst * z) .* dxOmegaCont(t, x);
    dXZKCont = cellfun(@(c) @(t,x,z) c * ones(size(x)), {0, 0; 0, 0}, 'UniformOutput', false);
    
    % Boundary conditions
    problemData.hDCont = problemData.hCont;
    problemData.gNCont = @(t,x,z) zeros(size(x));
    
    % Right hand side
    problemData.fCont = @(t,x,z) dThCont(t,x,z) - ...
                          dXZKCont{1,1}(t,x,z) .* dXhCont(t,x,z)  - problemData.KCont{1,1}(t,x,z) .* dXdXhCont(t,x,z) - ...
                          dXZKCont{1,2}(t,x,z) .* dZhCont(t,x,z)  - problemData.KCont{1,2}(t,x,z) .* dXdZhCont(t,x,z) - ...
                          dXZKCont{2,1}(t,x,z) .* dXhCont(t,x,z)  - problemData.KCont{2,1}(t,x,z) .* dXdZhCont(t,x,z) - ...
                          dXZKCont{2,2}(t,x,z) .* dZhCont(t,x,z)  - problemData.KCont{2,2}(t,x,z) .* dZdZhCont(t,x,z);
                   
  case 'convergence'
    % width and height of computational domain
    domainWidth = [0, 10];
    domainHeight = [0, 10];
    idDirichlet = [1, 2, 3, 4]; idNeumann = -1;
    % Analytical solution
    problemData.hCont = @(t,x,z) cos(x + t) .* cos(z + t);
    problemData.q1Cont = @(t,x,z) sin(x + t) .* cos(z + t);
    problemData.q2Cont = @(t,x,z) cos(x + t) .* sin(z + t);
    % Diffusion matrix
    problemData.KCont = { @(t,x,z) exp(z/5) , @(t,x,z) 0.5 * ones(size(x)) ; ...
                          @(t,x,z) 1/3 * ones(size(x)), @(t,x,z) exp(x/5) };
    % Derivatives
    dThCont = @(t,x,z) -sin(x + t) .* cos(z + t) - cos(x + t) .* sin(z + t);
    dXhCont = @(t,x,z) -sin(x + t) .* cos(z + t);
    dZhCont = @(t,x,z) -cos(x + t) .* sin(z + t);
    dXdXhCont = @(t,x,z) -cos(x + t) .* cos(z + t);
    dZdZhCont = @(t,x,z) -cos(x + t) .* cos(z + t);
    dXdZhCont = @(t,x,z) sin(x + t) .* sin(z + t);
    dXZKCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
                 @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
    % Boundary conditions
    problemData.hDCont = problemData.hCont;
    problemData.gNCont = @(t,x,z) zeros(size(x));
    % Right hand side
    problemData.fCont = @(t,x,z) dThCont(t,x,z) - ...
                          dXZKCont{1,1}(t,x,z) .* dXhCont(t,x,z)  - problemData.KCont{1,1}(t,x,z) .* dXdXhCont(t,x,z) - ...
                          dXZKCont{1,2}(t,x,z) .* dZhCont(t,x,z)  - problemData.KCont{1,2}(t,x,z) .* dXdZhCont(t,x,z) - ...
                          dXZKCont{2,1}(t,x,z) .* dXhCont(t,x,z)  - problemData.KCont{2,1}(t,x,z) .* dXdZhCont(t,x,z) - ...
                          dXZKCont{2,2}(t,x,z) .* dZhCont(t,x,z)  - problemData.KCont{2,2}(t,x,z) .* dZdZhCont(t,x,z);
end % switch

problemData = setdefault(problemData, 'isHotstart', isHotstart);
if problemData.isHotstart, problemData = setdefault(problemData, 'hotstartFile', hotstartFile); end

problemData = setdefault(problemData, 'idNeumann', idNeumann);
problemData = setdefault(problemData, 'idDirichlet', idDirichlet);

%% Domain and triangulation.
problemData.generateGrid = @(numElem) domainRectTrap(domainWidth, domainHeight, numElem);

% Boundary parts (0 = int, 1 = bot, 2 = right, 3 = top, 4 = left)
checkMultipleIds = @(idE0T, ids) logical(sum(bsxfun(@eq, idE0T, reshape(ids, 1, 1, length(ids))), 3));

problemData.generateMarkE0Tint = @(g) g.idE0T == 0;
problemData.generateMarkE0TbdrCoupling = @(g) g.idE0T == 3;
problemData.generateMarkE0TbdrN = @(g) checkMultipleIds(g.idE0T, problemData.idNeumann) & ~(problemData.isCoupling & g.idE0T == 3);
problemData.generateMarkE0TbdrD = @(g) checkMultipleIds(g.idE0T, problemData.idDirichlet) & ~(problemData.isCoupling & g.idE0T == 3);
end % function