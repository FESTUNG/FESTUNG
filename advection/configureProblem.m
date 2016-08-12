function problemData = configureProblem(problemData)
%% Parameters.
problemData.hmax        = 2^-5; % maximum edge length of triangle
problemData.p           = 0; % local polynomial degree
problemData.ordRK       = min(problemData.p+1,3); % order of Runge Kutta time stepper.
problemData.numSteps    = 320; % number of time steps
problemData.tEnd        = 1; % end time

problemData.isVisGrid   = false; % visualization of grid
problemData.isVisSol    = true; % visualization of solution
problemData.isSlopeLim  = false; % slope limiting
problemData.typeSlopeLim = ''; % Type of slope limiter (linear, hierarch_vert, strict)

problemData.outputFrequency = 32; % no visualization of every timestep
problemData.outputBasename  = ['output' filesep 'solution' problemData.typeSlopeLim]; % Basename of output files
problemData.outputTypes     = cellstr('vtk'); % solution output file types
%% Parameter check.
assert(problemData.p >= 0 && problemData.p <= 4, 'Polynomial order must be zero to four.')
assert(problemData.ordRK >= 1 && problemData.ordRK <= 3, 'Order of Runge Kutta must be zero to three.')
assert(problemData.hmax > 0, 'Maximum edge length must be positive.')
assert(problemData.numSteps > 0, 'Number of time steps must be positive.')
assert(~problemData.isSlopeLim || problemData.p > 0, 'Slope limiting only available for p > 0.')
%% Coefficients and boundary data
problemData.c0Cont = @(x1, x2)  sin(x1+x2)+1;
problemData.fCont = @(t,x1,x2) zeros(size(x1));
problemData.u1Cont = @(t,x1,x2) x1==x1;
problemData.u2Cont = @(t,x1,x2) 0*x1;
problemData.cDCont = @(t,x1,x2) sin(x1+x2-t)+1;
problemData.gNCont = @(t,x1,x2) zeros(size(x1));
end % function