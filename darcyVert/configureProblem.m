function problemData = configureProblem(problemData)

%% Parameters.
domainWidth = 100;  % width of computational domain
domainHeight = 2;  % height of computational domain
problemData.numElem = [2, 2];  % number of elements per direction
problemData.p = 1; % local polynomial degree
problemData.qOrd = 4; % order of quadrature formula
problemData.t0 = 0; % start time
problemData.tEnd = 0.1; % end time
problemData.numSteps = 20; % number of time steps
problemData.isVisGrid = true; % visualization of grid
problemData.isVisSol = true; % visualization of solution
problemData.eta = 1; % penalty parameter (eta>0)
problemData.outputBasename = ['output' filesep 'solution_darcy' ]; % Basename of output files
problemData.outputTypes = { 'vtk' };

%% Parameter check.
assert(problemData.p >= 0 && problemData.p <= 2, 'Polynomial order must be zero to two.')
assert(problemData.numSteps > 0, 'Number of time steps must be positive.')

%% Coefficients and boundary data.
S_0    = 1;
paramD = 0.01;
% Diffusion matrix
problemData.K11 = @(t,x1,x2) x1-x1 + exp(x2/5) / S_0;
problemData.K12 = @(t,x1,x2) x1-x1 + 0.5 / S_0;
problemData.K21 = @(t,x1,x2) x2-x2 + 1/3 / S_0;
problemData.K22 = @(t,x1,x2) x2-x2 + exp(x1/5) / S_0;
% Analytical solution
problemData.hCont = @(t,x1,x2) sin(paramD * (t+x1)) .* sin(paramD * (t+x2));
problemData.q1Cont = @(t,x1,x2) -paramD * cos(paramD * (t+x1)) .* sin(paramD * (t+x2));
problemData.q2Cont = @(t,x1,x2) -paramD * sin(paramD * (t+x1)) .* cos(paramD * (t+x2));
% Boundary conditions
problemData.gNCont = @(t,x1,x2) (2 * (x1 > 1) - 1) .* ( problemData.q1Cont(t,x1,x2) .* problemData.K11(t,x1,x2) ...
                       + problemData.q2Cont(t,x1,x2) .* problemData.K12(t,x1,x2) );        % Neumann boundary condition
problemData.hDCont = problemData.hCont; % Dirichlet boundary condition
% Analytical right hand side
problemData.f = @(t,x1,x2) ( ...
                paramD^2 * problemData.hCont(t,x1,x2) .* ( problemData.K11(t,x1,x2) + problemData.K22(t,x1,x2) ) - ...
                paramD^2 * cos(paramD * (t+x1)) .* cos(paramD * (t+x2)) + ...
                paramD * cos(paramD * (t+x1)) .* sin(paramD * (t+x2)) + ...
                paramD * sin(paramD * (t+x1)) .* cos(paramD * (t+x2)) ) / S_0;

%% Domain and triangulation.
problemData.generateGrid = @(numElem) domainRectTrap([0, domainWidth], [0, domainHeight], numElem);
end % function