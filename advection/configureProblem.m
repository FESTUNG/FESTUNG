function problemData = configureProblem(problemData)
%% Parameters.
problemData.hmax        = 2^-6; % maximum edge length of triangle
problemData.p           = 2; % local polynomial degree
problemData.ordRK       = min(problemData.p+1,3); % order of Runge Kutta time stepper.
problemData.numSteps    = 100; % number of time steps
problemData.tEnd        = (problemData.numSteps/3142)*2*pi; % end time

problemData.isVisGrid   = false; % visualization of grid
problemData.isVisSol    = true; % visualization of solution
problemData.isSlopeLim  = true; % slope limiting
problemData.typeSlopeLim = 'hierarch_vert'; % Type of slope limiter (linear, hierarch_vert, strict)

problemData.outputFrequency = 100; % no visualization of every timestep
problemData.outputBasename  = ['output' filesep 'solution_advection_' problemData.typeSlopeLim]; % Basename of output files
problemData.outputTypes     = {'vtk'}; % solution output file types
%% Parameter check.
assert(problemData.p >= 0 && problemData.p <= 4, 'Polynomial order must be zero to four.')
assert(problemData.ordRK >= 1 && problemData.ordRK <= 3, 'Order of Runge Kutta must be zero to three.')
assert(problemData.hmax > 0, 'Maximum edge length must be positive.')
assert(problemData.numSteps > 0, 'Number of time steps must be positive.')
assert(~problemData.isSlopeLim || problemData.p > 0, 'Slope limiting only available for p > 0.')
%% Coefficients and boundary data (LeVeque's solid body rotation).
G = @(x1, x2, x1_0, x2_0) (1/0.15) * sqrt((x1-x1_0).^2 + (x2-x2_0).^2);
problemData.c0Cont = @(x1, x2) ((x1 - 0.5).^2 + (x2 - 0.75).^2 <= 0.0225 & (x1 <= 0.475 | x1 >= 0.525 | x2 >= 0.85)) + ...
                    (1-G(x1, x2, 0.5, 0.25)) .* ((x1 - 0.5).^2 + (x2 - 0.25).^2 <= 0.0225) + ...
                    0.25*(1+cos(pi*G(x1, x2, 0.25, 0.5))).*((x1 - 0.25).^2 + (x2 - 0.5).^2 <= 0.0225);
problemData.fCont = @(t,x1,x2) zeros(size(x1));
problemData.u1Cont = @(t,x1,x2) 0.5 - x2;
problemData.u2Cont = @(t,x1,x2) x1 - 0.5;
problemData.cDCont = @(t,x1,x2) zeros(size(x1));
problemData.gNCont = @(t,x1,x2) zeros(size(x1));
%% Domain and triangulation configuration.
% Triangulate unit square using pdetool (if available or Friedrichs-Keller otherwise).
if license('checkout','PDE_Toolbox')
  problemData.generateGridData = @(hmax) domainPolygon([0 1 1 0], [0 0 1 1], hmax);
else
  fprintf('PDE_Toolbox not available. Using Friedrichs-Keller triangulation.\n');
  problemData.generateGridData = @domainSquare;
end % if
% Specify edge ids of boundary conditions
problemData.generateMarkE0Tint = @(g) g.idE0T == 0;
problemData.generateMarkE0TbdrN = @(g) false(g.numT,3);
problemData.generateMarkE0TbdrD = @(g) ~(g.markE0Tint | g.markE0TbdrN);
end % function