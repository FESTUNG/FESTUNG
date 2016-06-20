function problemData = configureProblem(problemData)

%% Parameters. 
problemData.p          = 2; % local polynomial degree
problemData.numSpecies = 2; % number of species involved

if ~isfield(problemData, 'hmax'), problemData.hmax        = 2^-6; end % maximum edge length of triangle
if ~isfield(problemData, 'ordRK'), problemData.ordRK       = min(problemData.p+1,3); end % order of Runge Kutta time stepper.
if ~isfield(problemData, 'tEnd'), problemData.tEnd        = 100/3142*2*pi; end % end time
if ~isfield(problemData, 'numSteps'), problemData.numSteps    = 100/3142*3142; end % number of time steps
if ~isfield(problemData, 'isVisGrid'), problemData.isVisGrid   = false; end % visualization of grid

%% Parameter check.
assert(problemData.p >= 0 && problemData.p <= 4, 'Polynomial order must be zero to four.')
assert(problemData.ordRK >= 1 && problemData.ordRK <= 3, 'Order of Runge Kutta must be zero to three.')
assert(problemData.hmax > 0, 'Maximum edge length must be positive.')
assert(problemData.numSteps > 0, 'Number of time steps must be positive.')

%% Coefficients and boundary data (LeVeque's solid body rotation).
G = @(x1, x2, x1_0, x2_0) (1/0.15) * sqrt((x1-x1_0).^2 + (x2-x2_0).^2);
if ~isfield(problemData, 'u1Cont'), problemData.u1Cont = @(t,x1,x2) 0.5 - x2; end % if
if ~isfield(problemData, 'u2Cont'), problemData.u2Cont = @(t,x1,x2) x1 - 0.5; end % if

problemData.isVisSol = cell(problemData.numSpecies,1);
problemData.isSlopeLim = cell(problemData.numSpecies,1);
problemData.typeSlopeLim = cell(problemData.numSpecies,1);
problemData.outputFrequency = cell(problemData.numSpecies,1);
problemData.outputBasename = cell(problemData.numSpecies,1);
problemData.outputTypes = cell(problemData.numSpecies,1);
problemData.c0Cont = cell(problemData.numSpecies,1);
problemData.fCont = cell(problemData.numSpecies,1);
problemData.cDCont = cell(problemData.numSpecies,1);
problemData.gNCont = cell(problemData.numSpecies,1);
for species = 1:problemData.numSpecies
  problemData.isVisSol{species}    = true; % visualization of solution
  problemData.isSlopeLim{species}  = true; % slope limiting
  problemData.typeSlopeLim{species} = 'hierarch_vert'; % Type of slope limiter (linear, hierarch_vert, strict)
  
  problemData.outputFrequency{species} = 100; % no visualization of every timestep
  problemData.outputBasename{species}  = ['solution_' num2str(species) '_' problemData.typeSlopeLim{species}]; % Basename of output files
  problemData.outputTypes{species}     = cellstr(['vtk';'tec']); % solution output file types
  
  %% Parameter check.
  assert(~problemData.isSlopeLim{species} || problemData.p > 0, 'Slope limiting only available for p > 0.')
  %% Coefficients and boundary data (LeVeque's solid body rotation).
  problemData.c0Cont{species} = @(x1, x2) (-1)^species*((x1 - 0.5).^2 + (x2 - 0.75).^2 <= 0.0225 & (x1 <= 0.475 | x1 >= 0.525 | x2 >= 0.85)) + ...
    (-1)^species*(1-G(x1, x2, 0.5, 0.25)) .* ((x1 - 0.5).^2 + (x2 - 0.25).^2 <= 0.0225) + ...
    (-1)^species*0.25*(1+cos(pi*G(x1, x2, 0.25, 0.5))).*((x1 - 0.25).^2 + (x2 - 0.5).^2 <= 0.0225);
  problemData.fCont{species} = @(t,x1,x2) zeros(size(x1));
  problemData.cDCont{species} = @(t,x1,x2) zeros(size(x1));
  problemData.gNCont{species} = @(t,x1,x2) zeros(size(x1));
end % for
end % function