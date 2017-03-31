function problemData = configureProblem(problemData)

%% Parameters.
% Name of testcase
problemData = setdefault(problemData, 'testcase', 'coupling');

problemData.eta = 1; % penalty parameter (eta>0)

% Number of elements in x- and y-direction
problemData = setdefault(problemData, 'numElem', [4, 2]);

% Local polynomial approximation order (0 to 5)
problemData = setdefault(problemData, 'p', 1);

% Order of quadrature rule
problemData = setdefault(problemData, 'qOrd', 2*problemData.p + 1);

% Time stepping parameters
problemData = setdefault(problemData, 't0', 0);  % start time
problemData = setdefault(problemData, 'tEnd', 0.1);  % end time
problemData = setdefault(problemData, 'numSteps', 10);  % number of time steps

% Discard time derivative and compute stationary solution
problemData = setdefault(problemData, 'isStationary', false);  

% Visualization settings
problemData = setdefault(problemData, 'isVisGrid', false);  % visualization of grid
problemData = setdefault(problemData, 'isVisSol', false);  % visualization of solution
problemData = setdefault(problemData, 'outputFrequency', 10); % no visualization of every timestep
problemData = setdefault(problemData, 'outputBasename', ...  % Basename of output files
                         ['output' filesep 'solution_darcy' ]); 
problemData = setdefault(problemData, 'outputTypes', { 'vtk' });  % Type of visualization files ('vtk, 'tec')

%% Parameter check.
assert(problemData.p >= 0 && problemData.p <= 5, 'Polynomial order must be zero to four.')
assert(problemData.numSteps > 0, 'Number of time steps must be positive.')

%% Coefficients and boundary data.
switch problemData.testcase
  case 'coupling'
    % width and height of computational domain
    domainWidth = [0, 100];
    domainHeight = [-2, 0];
    idBdrD = [1, 2, 4]; idBdrN = -1; idBdrCoupling = 3;
    % Analytical solution
    a = 0.1;
    b = 0.1;
    c = 0.5;
    d = 0.1;
    k = 1;
    problemData.hCont = @(t,x,z) a * cos(b*x + c*t) .* cos(d*z) + 3;
    problemData.q1Cont = @(t,x,z) a * b * sin(b*x + c*t) .* cos(d*z);
    problemData.q2Cont = @(t,x,z) a * d * cos(b*x + c*t) .* sin(d*z);
    % Diffusion matrix
    problemData.KCont = cellfun(@(c) @(t,x,z) c * ones(size(x)), {k, 0; 0, k}, 'UniformOutput', false);
    % Derivatives
    dThCont = @(t,x,z) -a * c * sin(b*x + c*t) .* cos(d*z);
    dXhCont = @(t,x,z) -a * b * sin(b*x + c*t) .* cos(d*z);
    dZhCont = @(t,x,z) -a * d * cos(b*x + c*t) .* sin(d*z);
    dXdXhCont = @(t,x,z) -a * b^2 * cos(b*x + c*t) .* cos(d*z);
    dZdZhCont = @(t,x,z) -a * d^2 * cos(b*x + c*t) .* cos(d*z);
    dXdZhCont = @(t,x,z) a * b * d * sin(b*x + c*t) .* sin(d*z);
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
                        
  case 'coupling2'
    % width and height of computational domain
    domainWidth = [0, 100];
    domainHeight = [0, 2];
    idBdrD = [1, 2, 4]; idBdrN = -1; idBdrCoupling = 3;
    % Analytical solution
    omega = 0.01;
    problemData.hCont = @(t,x,z) sin(omega * (t+x)) .* sin(omega * (t+z));
    problemData.q1Cont = @(t,x,z) -omega * cos(omega * (t+x)) .* sin(omega * (t+z));
    problemData.q2Cont = @(t,x,z) -omega * sin(omega * (t+x)) .* cos(omega * (t+z));
    % Diffusion matrix
    problemData.KCont = { @(t,x,z) exp(z/5) , @(t,x,z) 0.5 * ones(size(x)) ; ...
                          @(t,x,z) 1/3 * ones(size(x)), @(t,x,z) exp(x/5) };
    % Derivatives
    dThCont = @(t,x,z) omega * cos(omega * (t+x)) .* sin(omega * (t+z)) + omega * sin(omega * (t+x)) .* cos(omega * (t+z));
    dXhCont = @(t,x,z) omega * cos(omega * (t+x)) .* sin(omega * (t+z));
    dZhCont = @(t,x,z) omega * sin(omega * (t+x)) .* cos(omega * (t+z));
    dXdXhCont = @(t,x,z) -omega^2 * sin(omega * (t+x)) .* sin(omega * (t+z));
    dZdZhCont = @(t,x,z) -omega^2 * sin(omega * (t+x)) .* sin(omega * (t+z));
    dXdZhCont = @(t,x,z) omega * cos(omega * (t+x)) .* cos(omega * (t+z));
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
                        
  
  case 'convergence'
    % width and height of computational domain
    domainWidth = [0, 10];
    domainHeight = [0, 10];
    idBdrD = [1, 2, 3, 4]; idBdrN = -1; idBdrCoupling = -1;
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
    domainWidth = [0, 1];
    domainHeight = [0, 1];
    idBdrD = [1, 2, 3, 4]; idBdrN = -1; idBdrCoupling = -1;
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
end % switch

%% Domain and triangulation.
problemData.generateGrid = @(numElem) domainRectTrap(domainWidth, domainHeight, numElem);

% Boundary parts (0 = int, 1 = bot, 2 = right, 3 = top, 4 = left)
checkMultipleIds = @(idE0T, ids) logical(sum(bsxfun(@eq, idE0T, reshape(ids, 1, 1, length(ids))), 3));

problemData.generateMarkE0Tint = @(g) g.idE0T == 0;
problemData.generateMarkE0TbdrN = @(g) checkMultipleIds(g.idE0T, idBdrN);
problemData.generateMarkE0TbdrD = @(g) checkMultipleIds(g.idE0T, idBdrD);
problemData.generateMarkE0TbdrCoupling = @(g) checkMultipleIds(g.idE0T, idBdrCoupling);
end % function