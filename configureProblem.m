function problemData = configureProblem(problemData)

%% Parameters.
problemData.width           = 100;
problemData.heightDarcy     = 2;
problemData.NX              = 16;         % maximum edge length of triangle
problemData.NZdarcy         = 16;         % 
problemData.p               = 1;          % local polynomial degree
problemData.ord             = 4;          % order of quadrature formula
problemData.tEnd            = 0.1;         % end time
problemData.numSteps        = 20;         % number of time steps
problemData.isVisSol        = false;      % visualization of solution
problemData.eta             = 1;          % penalty parameter (eta>0)
problemData.outputBasename  = 'solution'; % Basename of output files
%problemData.outputTypes     = cellstr(['vtk']);

%% Parameter check.
assert(problemData.p >= 0 && problemData.p <= 2, 'Polynomial order must be zero to two.')
assert(problemData.numSteps > 0, 'Number of time steps must be positive.')

%% Coefficients and boundary data.
problemData.S_0     = 1;
problemData.paramD  = 0.1;
problemData.K11     = @(t,x,y) x-x + exp(y/5) / S_0;
problemData.K12     = @(t,x,y) x-x + 1/2 / S_0;
problemData.K21     = @(t,x,y) y-y + 1/3 / S_0;
problemData.K22     = @(t,x,y) y-y + exp(x/5) / S_0;
problemData.h_D     = @(t,x,y) sin(parD * (t+x)) .* sin(parD * (t+y));     % Dirichlet boundary condition
problemData.q1      = @(t,x,y) -parD * cos(parD * (t+x)) .* sin(parD * (t+y));
problemData.q2      = @(t,x,y) -parD * sin(parD * (t+x)) .* cos(parD * (t+y));
problemData.g_N     = @(t,x,y) (2 * (x > 1) - 1) .* ( q1SolTimePM(t,x,y) .* K11TimePM(t,x,y) ...
                        + q2SolTimePM(t,x,y) .* K12TimePM(t,x,y) );        % Neumann boundary condition
problemData.hSol    = problemdata.h_D;
problemData.f       = @(t,x,y) ( parD^2 * hSolTimePM(t,x,y) .* ( K11TimePM(t,x,y) + K22TimePM(t,x,y) ) ...
                        - parD^2 * cos(parD * (t+x)) .* cos(parD * (t+y)) ...
                        + ( parD * cos(parD * (t+x)) .* sin(parD * (t+y)) + parD * sin(parD * (t+x)) .* cos(parD * (t+y)) ) ) / S_0;

end % function