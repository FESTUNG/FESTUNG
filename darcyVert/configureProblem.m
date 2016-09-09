function problemData = configureProblem(problemData)

%% Parameters.
domainWidth = 100;  % width of computational domain
domainHeight = 2;  % height of computational domain
problemData.numElem = [16, 16];  % number of elements per direction
problemData.p = 1; % local polynomial degree
problemData.qOrd = 4; % order of quadrature formula
problemData.t0 = 0; % start time
problemData.tEnd = 0.1; % end time
problemData.numSteps = 20; % number of time steps
problemData.isVisGrid = true; % visualization of grid
problemData.isVisSol = false; % visualization of solution
problemData.eta = 1; % penalty parameter (eta>0)
problemData.outputBasename = ['output' filesep 'solution_darcy' ]; % Basename of output files
problemData.outputTypes = { 'vtk' };

%% Parameter check.
assert(problemData.p >= 0 && problemData.p <= 2, 'Polynomial order must be zero to two.')
assert(problemData.numSteps > 0, 'Number of time steps must be positive.')

%% Coefficients and boundary data.
S_0    = 1;
paramD = 0.1;
problemData.K11    = @(t,x,y) x-x + exp(y/5) / S_0;
problemData.K12    = @(t,x,y) x-x + 1/2 / S_0;
problemData.K21    = @(t,x,y) y-y + 1/3 / S_0;
problemData.K22    = @(t,x,y) y-y + exp(x/5) / S_0;
problemData.h_D    = @(t,x,y) sin(paramD * (t+x)) .* sin(paramD * (t+y));     % Dirichlet boundary condition
problemData.q1     = @(t,x,y) -paramD * cos(paramD * (t+x)) .* sin(paramD * (t+y));
problemData.q2     = @(t,x,y) -paramD * sin(paramD * (t+x)) .* cos(paramD * (t+y));
problemData.g_N    = @(t,x,y) (2 * (x > 1) - 1) .* ( q1SolTimePM(t,x,y) .* K11TimePM(t,x,y) ...
                       + q2SolTimePM(t,x,y) .* K12TimePM(t,x,y) );        % Neumann boundary condition
problemData.hSol   = domainHeight;
problemData.f      = @(t,x,y) ( paramD^2 * hSolTimePM(t,x,y) .* ( K11TimePM(t,x,y) + K22TimePM(t,x,y) ) ...
                       - paramD^2 * cos(paramD * (t+x)) .* cos(paramD * (t+y)) ...
                       + ( paramD * cos(paramD * (t+x)) .* sin(paramD * (t+y)) ...
                       + paramD * sin(paramD * (t+x)) .* cos(paramD * (t+y)) ) ) / S_0;

%% Domain
problemData.generateGrid = @(nx, nz) domainRectTrap([0, domainWidth], [0, domainHeight], [nx, nz]);
end % function