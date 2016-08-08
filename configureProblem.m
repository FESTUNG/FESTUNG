function problemData = configureProblem(problemData)
%% Parameters. 
% Set default values if they are not yet available in problemData
problemData = setdefault(problemData, 'numSpecies', 1);  % number of transported species
problemData = setdefault(problemData, 'p'         , 2);  % local polynomial degree (TODO: allow different approximation orders for each species)
problemData = setdefault(problemData, 'hmax'      , 2^-6);  % maximum edge length of triangle
problemData = setdefault(problemData, 'ordRK'     , min(problemData.p+1,3));  % order of Runge Kutta time stepper
problemData = setdefault(problemData, 'numSteps'  , 3142);  % number of time steps
problemData = setdefault(problemData, 'tEnd'      , (problemData.numSteps/3142)*2*pi);  % end time
problemData = setdefault(problemData, 'isVisGrid' , false);  % visualization of grid

%% Parameter check.
assert(problemData.p >= 0 && problemData.p <= 4, 'Polynomial order must be zero to four.')
assert(problemData.ordRK >= 1 && problemData.ordRK <= 3, 'Order of Runge Kutta must be zero to three.')
assert(problemData.hmax > 0, 'Maximum edge length must be positive.')
assert(problemData.numSteps > 0, 'Number of time steps must be positive.')

%% Coefficients and boundary data (LeVeque's solid body rotation).
G = @(x1, x2, x1_0, x2_0) (1/0.15) * sqrt((x1-x1_0).^2 + (x2-x2_0).^2);
c0Cont = @(x1, x2) ((x1 - 0.5).^2 + (x2 - 0.75).^2 <= 0.0225 & (x1 <= 0.475 | x1 >= 0.525 | x2 >= 0.85)) + ...
                    (1-G(x1, x2, 0.5, 0.25)) .* ((x1 - 0.5).^2 + (x2 - 0.25).^2 <= 0.0225) + ...
                    0.25*(1+cos(pi*G(x1, x2, 0.25, 0.5))).*((x1 - 0.25).^2 + (x2 - 0.5).^2 <= 0.0225);
problemData = setdefault(problemData, 'u1Cont', @(t,x1,x2) 0.5 - x2);
problemData = setdefault(problemData, 'u2Cont', @(t,x1,x2) x1 - 0.5);

problemData = setdefault(problemData, 'outputBasename', ['output' filesep 'solution_transport']); % Basename of output files
problemData = setdefault(problemData, 'outputTypes', {'vtk'}); % solution output file types

problemData.isVisSol = cell(problemData.numSpecies,1);
problemData.isSlopeLim = cell(problemData.numSpecies,1);
problemData.typeSlopeLim = cell(problemData.numSpecies,1);
problemData.outputFrequency = cell(problemData.numSpecies,1);
problemData.c0Cont = cell(problemData.numSpecies,1);
problemData.fCont = cell(problemData.numSpecies,1);
problemData.cDCont = cell(problemData.numSpecies,1);
problemData.gNCont = cell(problemData.numSpecies,1);
problemData.reactions = cell(problemData.numSpecies,1);

for species = 1:problemData.numSpecies
  problemData.isVisSol{species}    = true; % visualization of solution
  problemData.isSlopeLim{species}  = true; % slope limiting
  problemData.typeSlopeLim{species} = 'hierarch_vert'; % Type of slope limiter (linear, hierarch_vert, strict)
  problemData.outputFrequency{species} = 100; % no visualization of every timestep
  
  %% Parameter check.
  assert(~problemData.isSlopeLim{species} || problemData.p > 0, 'Slope limiting only available for p > 0.')
  
  %% Coefficients and boundary data (LeVeque's solid body rotation).
  problemData.c0Cont{species} = @(x1, x2) species * c0Cont(x1, x2);
  problemData.fCont{species} = @(t,x1,x2) zeros(size(x1));
  problemData.cDCont{species} = @(t,x1,x2) zeros(size(x1));
  problemData.gNCont{species} = @(t,x1,x2) zeros(size(x1));
end % for

%% Reaction term definitions.
% NPZ model
% parameters (Notation as in paper)
I0 = 1;
Vm = 1;
ks = 1;
Rm = 1;
ep = 1;
gamma = 1;
I = @(t,x1,x2) (x1==x1);
f = @(t,x1,x2) I(t,x1,x2) /I0; % linear response
g = @(t,x1,x2,N) Vm ./ (ks + N); % Michaelis-Menten uptake
h = @(t,x1,x2,P) Rm * P; % linear grazing
i = @(t,x1,x2,P) ep; % linear death rate
j = @(t,x1,x2,Z) ep; % linear death rate

problemData.reactions{1} = @(t,x1,x2,c) 0*x1;
problemData.reactions{2} = @(t,x1,x2,c) 0*x1;
problemData.reactions{3} = @(t,x1,x2,c) 0*x1;
% problemData.reactions{1} = @(t,x1,x2,c) f(t,x1,x2) .* g(t,x1,x2,c{3}) .* c{1} - h(t,x1,x2,c{1}) .* c{2} - i(t,x1,x2,c{1}) .* c{1};
% problemData.reactions{2} = @(t,x1,x2,c) gamma * h(t,x1,x2,c{1}) .* c{2} - j(t,x1,x2,c{2}) .* c{2};
% problemData.reactions{3} = @(t,x1,x2,c) -f(t,x1,x2) .* g(t,x1,x2,c{3}) .* c{1} + (1-gamma) .* h(t,x1,x2,c{1}) .* c{2} + i(t,x1,x2,c{1}) .* c{1} + j(t,x1,x2,c{2}) .* c{2};
end % function