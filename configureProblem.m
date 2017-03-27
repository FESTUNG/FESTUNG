function problemData = configureProblem(problemData)

%% Parameters.
% Name of testcase
problemData = setdefault(problemData, 'testcase', 'convergence');

problemData.eta = 1; % penalty parameter (eta>0)

% Number of elements in x- and y-direction
problemData = setdefault(problemData, 'numElem', [32, 16]);

% Local polynomial approximation order (0 to 5)
problemData = setdefault(problemData, 'p', 1);

% Order of quadrature rule
problemData = setdefault(problemData, 'qOrd', 2*problemData.p + 1);

% Time stepping parameters
problemData = setdefault(problemData, 't0', 0);  % start time
problemData = setdefault(problemData, 'tEnd', 1);  % end time
problemData = setdefault(problemData, 'numSteps', 40);  % number of time steps

% Discard time derivative and compute stationary solution
problemData = setdefault(problemData, 'isStationary', false);  

% Visualization settings
problemData = setdefault(problemData, 'isVisGrid', false);  % visualization of grid
problemData = setdefault(problemData, 'isVisSol', true);  % visualization of solution
problemData = setdefault(problemData, 'outputFrequency', 10); % no visualization of every timestep
problemData = setdefault(problemData, 'outputBasename', ...  % Basename of output files
                         ['output' filesep 'solution_darcy' ]); 
problemData = setdefault(problemData, 'outputTypes', { 'vtk' });  % Type of visualization files ('vtk, 'tec')

%% Parameter check.
assert(problemData.p >= 0 && problemData.p <= 5, 'Polynomial order must be zero to four.')
assert(problemData.numSteps > 0, 'Number of time steps must be positive.')

%% Coefficients and boundary data.
switch problemData.testcase
  case 'convergence'
    % width and height of computational domain
    domainWidth = 10;
    domainHeight = 10;
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
                        
  case 'convergence2'
    % width and height of computational domain
    domainWidth = 1;
    domainHeight = 1;
    % Analytical solution
    problemData.hCont = @(t,x,z) cos(7 * x) .* cos(7 * z);
    problemData.q1Cont = @(t,x,z) 7 * sin(x) .* cos(z);
    problemData.q2Cont = @(t,x,z) 7 * cos(x) .* sin(z);
    % Diffusion matrix
    problemData.KCont = { @(t,x,z) exp(x+z) , @(t,x,z) zeros(size(x)) ; ...
                          @(t,x,z) zeros(size(x)), @(t,x,z) exp(x+z) };
    % Derivatives
    dThCont = @(t,x,z) zeros(size(x));
    dXhCont = @(t,x,z) -7 * sin(7 * x) .* cos(7 * z);
    dZhCont = @(t,x,z) -7 * cos(7 * x) .* sin(7 * z);
    dXdXhCont = @(t,x,z) -49 * cos(7 * x) .* cos(7 * z);
    dZdZhCont = @(t,x,z) -49 * cos(7 * x) .* cos(7 * z);
    dXdZhCont = @(t,x,z) 49 * sin(7 * x) .* sin(7 * z);
    dXZKCont = problemData.KCont;
    % Boundary conditions
    problemData.hDCont = problemData.hCont;
    problemData.gNCont = @(t,x,z) zeros(size(x));
    % Right hand side
    problemData.fCont = @(t,x,z) dThCont(t,x,z) - ...
                          dXZKCont{1,1}(t,x,z) .* dXhCont(t,x,z)  - problemData.KCont{1,1}(t,x,z) .* dXdXhCont(t,x,z) - ...
                          dXZKCont{1,2}(t,x,z) .* dZhCont(t,x,z)  - problemData.KCont{1,2}(t,x,z) .* dXdZhCont(t,x,z) - ...
                          dXZKCont{2,1}(t,x,z) .* dXhCont(t,x,z)  - problemData.KCont{2,1}(t,x,z) .* dXdZhCont(t,x,z) - ...
                          dXZKCont{2,2}(t,x,z) .* dZhCont(t,x,z)  - problemData.KCont{2,2}(t,x,z) .* dZdZhCont(t,x,z);

  case 'coupling'
    % width and height of computational domain
    domainWidth = 10;
    domainHeight = 10;
    S_0    = 1;
    paramD = 0.01;
    % Diffusion matrix
    problemData.KCont = cell(2,2);
    problemData.KCont{1,1} = @(t,x1,x2) x1-x1 + exp(x2/5) / S_0;
    problemData.KCont{1,2} = @(t,x1,x2) x1-x1 + 0.5 / S_0;
    problemData.KCont{2,1} = @(t,x1,x2) x2-x2 + 1/3 / S_0;
    problemData.KCont{2,2} = @(t,x1,x2) x2-x2 + exp(x1/5) / S_0;
    % Analytical solution
    problemData.hCont = @(t,x1,x2) sin(paramD * (t+x1)) .* sin(paramD * (t+x2));
    problemData.q1Cont = @(t,x1,x2) -paramD * cos(paramD * (t+x1)) .* sin(paramD * (t+x2));
    problemData.q2Cont = @(t,x1,x2) -paramD * sin(paramD * (t+x1)) .* cos(paramD * (t+x2));
    % Boundary conditions
    problemData.gNCont = @(t,x1,x2) (2 * (x1 > 1) - 1) .* ( problemData.q1Cont(t,x1,x2) .* problemData.KCont{1,1}(t,x1,x2) ...
                           + problemData.q2Cont(t,x1,x2) .* problemData.KCont{1,2}(t,x1,x2) );        % Neumann boundary condition
    problemData.hDCont = problemData.hCont; % Dirichlet boundary condition
    % Analytical right hand side
    problemData.fCont = @(t,x1,x2) ( ...
                        paramD^2 * problemData.hCont(t,x1,x2) .* ( problemData.KCont{1,1}(t,x1,x2) + problemData.KCont{2,2}(t,x1,x2) ) - ...
                        paramD^2 * cos(paramD * (t+x1)) .* cos(paramD * (t+x2)) + ...
                        paramD * cos(paramD * (t+x1)) .* sin(paramD * (t+x2)) + ...
                        paramD * sin(paramD * (t+x1)) .* cos(paramD * (t+x2)) ) / S_0;
end % switch

%% Domain and triangulation.
problemData.generateGrid = @(numElem) domainRectTrap([0, domainWidth], [0, domainHeight], numElem);

% Boundary parts (0 = int, 1 = bot, 2 = right, 3 = top, 4 = left)
checkMultipleIds = @(idE0T, ids) logical(sum(bsxfun(@eq, idE0T, reshape(ids, 1, 1, length(ids))), 3));

problemData.generateMarkE0Tint = @(g) g.idE0T == 0;
problemData.generateMarkE0TbdrN = @(g) checkMultipleIds(g.idE0T, -1);
problemData.generateMarkE0TbdrD = @(g) checkMultipleIds(g.idE0T, [1 2 3 4]);
end % function