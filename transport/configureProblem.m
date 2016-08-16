function problemData = configureProblem(problemData)

%% Configuration to use: 
% - 'rotation' calls configureRotation()
% - 'analytical' calls configureAnalyticalTest()
problemData = setdefault(problemData, 'configSource', 'analytical');

%% Parameters. 
% Set default values if they are not yet available in problemData
problemData = setdefault(problemData, 'p'         , 0);  % local polynomial degree (TODO: allow different approximation orders for each species)
problemData = setdefault(problemData, 'ordRK'     , min(problemData.p+1,3));  % order of Runge Kutta time stepper
problemData = setdefault(problemData, 'isVisGrid' , false);  % visualization of grid

%% Parameter check.
assert(problemData.p >= 0 && problemData.p <= 4, 'Polynomial order must be zero to four.')
assert(problemData.ordRK >= 1 && problemData.ordRK <= 3, 'Order of Runge Kutta must be zero to three.')

problemData = setdefault(problemData, 'numSpecies', 3);  % number of transported species

problemData.isVisSol = cell(problemData.numSpecies,1);
problemData.isSlopeLim = cell(problemData.numSpecies,1);
problemData.typeSlopeLim = cell(problemData.numSpecies,1);
problemData.outputFrequency = cell(problemData.numSpecies,1);
problemData.outputBasename = cell(problemData.numSpecies,1);
problemData.outputTypes = cell(problemData.numSpecies,1);
problemData.cH0Cont = cell(problemData.numSpecies,1);
problemData.fCont = cell(problemData.numSpecies,1);
problemData.cDCont = cell(problemData.numSpecies,1);
problemData.gNCont = cell(problemData.numSpecies,1);
problemData.reactions = cell(problemData.numSpecies,1);

%% Simulation scenario specific parameters
switch problemData.configSource
  case 'rotation'
    problemData.isSolutionAvailable = false;
    problemData = setdefault(problemData, 'hmax'      , 2^-6);  % maximum edge length of triangle
    problemData = setdefault(problemData, 'numSteps'  , 100);  % number of time steps
    problemData = setdefault(problemData, 'tEnd'      , (problemData.numSteps/3142)*2*pi);  % end time
    problemData = configureRotation(problemData);
  case 'analytical'
    problemData.isSolutionAvailable = true;
    problemData = setdefault(problemData, 'hmax'      , 1);  % maximum edge length of triangle
    problemData = setdefault(problemData, 'numSteps'  , 100);  % number of time steps
    problemData = setdefault(problemData, 'tEnd'      , 1);  % end time
    problemData = configureAnalyticalTest(problemData);
  otherwise
    error('Invalid config source.')
end % switch

assert(problemData.hmax > 0, 'Maximum edge length must be positive.')
assert(problemData.numSteps > 0, 'Number of time steps must be positive.')

end % function

%% LeVeque's solid body rotation
function problemData = configureRotation(problemData)

problemData = setdefault(problemData, 'hCont', @(t,x1,x2) x1==x1); % TODO order of arguments
problemData = setdefault(problemData, 'uCont', @(t,x1,x2) 0.5 - x2);
problemData = setdefault(problemData, 'vCont', @(t,x1,x2) x1 - 0.5);
problemData = setdefault(problemData, 'uHCont', @(t,x1,x2) problemData.hCont(t,x1,x2) .* problemData.uCont(t,x1,x2));
problemData = setdefault(problemData, 'vHCont', @(t,x1,x2) problemData.hCont(t,x1,x2) .* problemData.vCont(t,x1,x2));

%% Coefficients and boundary data (LeVeque's solid body rotation).
G = @(x1, x2, x1_0, x2_0) (1/0.15) * sqrt((x1-x1_0).^2 + (x2-x2_0).^2);
cH0Cont = @(x1, x2) ((x1 - 0.5).^2 + (x2 - 0.75).^2 <= 0.0225 & (x1 <= 0.475 | x1 >= 0.525 | x2 >= 0.85)) + ...
                    (1-G(x1, x2, 0.5, 0.25)) .* ((x1 - 0.5).^2 + (x2 - 0.25).^2 <= 0.0225) + ...
                    0.25*(1+cos(pi*G(x1, x2, 0.25, 0.5))).*((x1 - 0.25).^2 + (x2 - 0.5).^2 <= 0.0225);
                  
for species = 1:problemData.numSpecies
  problemData.isVisSol{species}    = true; % visualization of solution
  problemData.isSlopeLim{species}  = false; % slope limiting
  problemData.typeSlopeLim{species} = 'hierarch_vert'; % Type of slope limiter (linear, hierarch_vert, strict)
  problemData.outputFrequency{species} = 100; % no visualization of every timestep
  
  %% Parameter check.
  assert(~problemData.isSlopeLim{species} || problemData.p > 0, 'Slope limiting only available for p > 0.')

  problemData.cH0Cont{species} = @(x1,x2) species * cH0Cont(x1, x2);
  problemData.cDCont{species} = @(t,x1,x2) zeros(size(x1));
  problemData.gNCont{species} = @(t,x1,x2) zeros(size(x1));
end % for

problemData.reactions{1} = @(t,x1,x2,c,cH) 0*x1;
problemData.reactions{2} = @(t,x1,x2,c,cH) 0*x1;
problemData.reactions{3} = @(t,x1,x2,c,cH) 0*x1;

problemData.fCont{1} = @(t,x1,x2) 0*x1;
problemData.fCont{2} = @(t,x1,x2) 0*x1;
problemData.fCont{3} = @(t,x1,x2) 0*x1;

end % function

%% Analytical solution
function problemData = configureAnalyticalTest(problemData)

A = 0.00001;
B = 0.00001;
C = 0.00001;
height = 0.00025;

problemData.zb_xCont = @(x1,x2) height*(x1==x1);
problemData.zb_yCont = @(x1,x2) height*(x1==x1);

problemData.hCont = @(t,x1,x2) C*(sin(0.0001*(x1-t)) + sin(0.0001*(x2-t))) - height*(x1+x2) + 10;
problemData.h_tCont = @(t,x1,x2) -0.0001*C*(cos(0.0001*(x1-t)) + cos(0.0001*(x2-t)));
problemData.h_xCont = @(t,x1,x2) 0.0001*C*cos(0.0001*(x1-t)) - problemData.zb_xCont(x1,x2);
problemData.h_yCont = @(t,x1,x2) 0.0001*C*cos(0.0001*(x2-t)) - problemData.zb_yCont(x1,x2);

problemData.uCont = @(t,x1,x2) A*sin(0.0001*(x1-t));
problemData.u_xCont = @(t,x1,x2) 0.0001*A*cos(0.0001*(x1-t));
problemData.vCont = @(t,x1,x2) B*sin(0.0001*(x2-t));
problemData.v_yCont = @(t,x1,x2) 0.0001*B*cos(0.0001*(x2-t));

problemData = setdefault(problemData, 'uHCont', @(t,x1,x2) problemData.hCont(t,x1,x2) .* problemData.uCont(t,x1,x2));
problemData = setdefault(problemData, 'vHCont', @(t,x1,x2) problemData.hCont(t,x1,x2) .* problemData.vCont(t,x1,x2));

% analytical solution
problemData.solCont = {@(t,x1,x2) sin(0.0001*(x1+x2-t))+2, @(t,x1,x2) cos(0.0001*(x1+x2-t))+2, @(t,x1,x2) -cos(0.0001*(x1+x2-t))+2};
sol_tCont = {@(t,x1,x2) -0.0001*cos(0.0001*(x1+x2-t)), @(t,x1,x2) 0.0001*sin(0.0001*(x1+x2-t)), @(t,x1,x2) -0.0001*sin(0.0001*(x1+x2-t))};
sol_xCont = {@(t,x1,x2) 0.0001*cos(0.0001*(x1+x2-t)), @(t,x1,x2) -0.0001*sin(0.0001*(x1+x2-t)), @(t,x1,x2) 0.0001*sin(0.0001*(x1+x2-t))};
sol_yCont = {@(t,x1,x2) 0.0001*cos(0.0001*(x1+x2-t)), @(t,x1,x2) -0.0001*sin(0.0001*(x1+x2-t)), @(t,x1,x2) 0.0001*sin(0.0001*(x1+x2-t))};

for species = 1:problemData.numSpecies
  problemData.isVisSol{species}    = true; % visualization of solution
  problemData.isSlopeLim{species}  = false; % slope limiting
  problemData.typeSlopeLim{species} = 'hierarch_vert'; % Type of slope limiter (linear, hierarch_vert, strict)
  
  problemData.outputFrequency{species} = problemData.numSteps/10; % no visualization of every timestep
  problemData.outputBasename{species}  = ['output' filesep 'p=' num2str(problemData.p) '_solution_' num2str(species)]; % Basename of output files
  problemData.outputTypes{species}     = cellstr('vtk'); % solution output file types
  
  %% Parameter check.
  assert(~problemData.isSlopeLim{species} || problemData.p > 0, 'Slope limiting only available for p > 0.')
  
  problemData.cH0Cont{species} = @(x1,x2) problemData.solCont{species}(0,x1,x2) .* problemData.hCont(0,x1,x2);
  problemData.cDCont{species} = @(t,x1,x2) problemData.solCont{species}(t,x1,x2);
  problemData.gNCont{species} = @(t,x1,x2) 0*x1; % TODO allow to neglect this if no Neumann boundary is used
end % for
  
%% Reaction term definitions.
% NPZ model
% parameters (Notation as in paper)
I0 = 1;
Vm = 0.0000001;
ks = 1;
Rm = 0.0000001;
ep = 0.0000001;
gamma = 0.5;
I = @(t,x1,x2) (x1==x1);
f = @(t,x1,x2) I(t,x1,x2) /I0; % linear response
g = @(t,x1,x2,N) Vm ./ (ks + N); % Michaelis-Menten uptake
h = @(t,x1,x2,P) Rm * P; % linear grazing
i = @(t,x1,x2,P) ep; % linear death rate
j = @(t,x1,x2,Z) ep; % linear death rate

problemData.reactions{1} = @(t,x1,x2,c,cH) f(t,x1,x2) .* g(t,x1,x2,c{3}) .* cH{1} - h(t,x1,x2,c{1}) .* cH{2} - i(t,x1,x2,c{1}) .* cH{1};
problemData.reactions{2} = @(t,x1,x2,c,cH) gamma * h(t,x1,x2,c{1}) .* cH{2} - j(t,x1,x2,c{2}) .* cH{2};
problemData.reactions{3} = @(t,x1,x2,c,cH) -f(t,x1,x2) .* g(t,x1,x2,c{3}) .* cH{1} + (1-gamma) .* h(t,x1,x2,c{1}) .* cH{2} + i(t,x1,x2,c{1}) .* cH{1} + j(t,x1,x2,c{2}) .* cH{2};

%% Attention: due to inconsistency the arguments t,x1,x2 appear in different order in swe and transport TODO
problemData.fCont{1} = @(t,x1,x2) sol_tCont{1}(t,x1,x2) .* problemData.hCont(t,x1,x2) + problemData.solCont{1}(t,x1,x2) .* problemData.h_tCont(t,x1,x2) + ... 
  problemData.u_xCont(t,x1,x2) .* problemData.hCont(t,x1,x2) .* problemData.solCont{1}(t,x1,x2) + problemData.uCont(t,x1,x2) .* problemData.h_xCont(t,x1,x2) .* problemData.solCont{1}(t,x1,x2) + problemData.uCont(t,x1,x2) .* problemData.hCont(t,x1,x2) .* sol_xCont{1}(t,x1,x2) + ... 
  problemData.v_yCont(t,x1,x2) .* problemData.hCont(t,x1,x2) .* problemData.solCont{1}(t,x1,x2) + problemData.vCont(t,x1,x2) .* problemData.h_yCont(t,x1,x2) .* problemData.solCont{1}(t,x1,x2) + problemData.vCont(t,x1,x2) .* problemData.hCont(t,x1,x2) .* sol_yCont{1}(t,x1,x2) - ... 
  (f(t,x1,x2) .* g(t,x1,x2,problemData.solCont{3}(t,x1,x2)) .* problemData.solCont{1}(t,x1,x2) - h(t,x1,x2,problemData.solCont{1}(t,x1,x2)) .* problemData.solCont{2}(t,x1,x2) - i(t,x1,x2,problemData.solCont{1}(t,x1,x2)) .* problemData.solCont{1}(t,x1,x2)) .* problemData.hCont(t,x1,x2);
problemData.fCont{2} = @(t,x1,x2) sol_tCont{2}(t,x1,x2) .* problemData.hCont(t,x1,x2) + problemData.solCont{2}(t,x1,x2) .* problemData.h_tCont(t,x1,x2) + ... 
  problemData.u_xCont(t,x1,x2) .* problemData.hCont(t,x1,x2) .* problemData.solCont{2}(t,x1,x2) + problemData.uCont(t,x1,x2) .* problemData.h_xCont(t,x1,x2) .* problemData.solCont{2}(t,x1,x2) + problemData.uCont(t,x1,x2) .* problemData.hCont(t,x1,x2) .* sol_xCont{2}(t,x1,x2) + ... 
  problemData.v_yCont(t,x1,x2) .* problemData.hCont(t,x1,x2) .* problemData.solCont{2}(t,x1,x2) + problemData.vCont(t,x1,x2) .* problemData.h_yCont(t,x1,x2) .* problemData.solCont{2}(t,x1,x2) + problemData.vCont(t,x1,x2) .* problemData.hCont(t,x1,x2) .* sol_yCont{2}(t,x1,x2) - ... 
  (gamma * h(t,x1,x2,problemData.solCont{1}(t,x1,x2)) .* problemData.solCont{2}(t,x1,x2) - j(t,x1,x2,problemData.solCont{2}(t,x1,x2)) .* problemData.solCont{2}(t,x1,x2)) .* problemData.hCont(t,x1,x2);
problemData.fCont{3} = @(t,x1,x2) sol_tCont{3}(t,x1,x2) .* problemData.hCont(t,x1,x2) + problemData.solCont{3}(t,x1,x2) .* problemData.h_tCont(t,x1,x2) + ... 
  problemData.u_xCont(t,x1,x2) .* problemData.hCont(t,x1,x2) .* problemData.solCont{3}(t,x1,x2) + problemData.uCont(t,x1,x2) .* problemData.h_xCont(t,x1,x2) .* problemData.solCont{3}(t,x1,x2) + problemData.uCont(t,x1,x2) .* problemData.hCont(t,x1,x2) .* sol_xCont{3}(t,x1,x2) + ... 
  problemData.v_yCont(t,x1,x2) .* problemData.hCont(t,x1,x2) .* problemData.solCont{3}(t,x1,x2) + problemData.vCont(t,x1,x2) .* problemData.h_yCont(t,x1,x2) .* problemData.solCont{3}(t,x1,x2) + problemData.vCont(t,x1,x2) .* problemData.hCont(t,x1,x2) .* sol_yCont{3}(t,x1,x2) - ... 
  (-f(t,x1,x2) .* g(t,x1,x2,problemData.solCont{3}(t,x1,x2)) .* problemData.solCont{1}(t,x1,x2) + (1-gamma) .* h(t,x1,x2,problemData.solCont{1}(t,x1,x2)) .* problemData.solCont{2}(t,x1,x2) + i(t,x1,x2,problemData.solCont{1}(t,x1,x2)) .* problemData.solCont{1}(t,x1,x2) + j(t,x1,x2,problemData.solCont{2}(t,x1,x2)) .* problemData.solCont{2}(t,x1,x2)) .* problemData.hCont(t,x1,x2);

% problemData.reactions{1} = @(t,x1,x2,c,cH) 0*x1;
% problemData.reactions{2} = @(t,x1,x2,c,cH) 0*x1;
% problemData.reactions{3} = @(t,x1,x2,c,cH) 0*x1;
 
% problemData.fCont{1} = @(t,x1,x2) sol_tCont{1}(t,x1,x2) .* problemData.hCont(t,x1,x2) + problemData.solCont{1}(t,x1,x2) .* problemData.h_tCont(t,x1,x2);
% problemData.fCont{2} = @(t,x1,x2) sol_tCont{2}(t,x1,x2) .* problemData.hCont(t,x1,x2) + problemData.solCont{2}(t,x1,x2) .* problemData.h_tCont(t,x1,x2);
% problemData.fCont{3} = @(t,x1,x2) sol_tCont{3}(t,x1,x2) .* problemData.hCont(t,x1,x2) + problemData.solCont{3}(t,x1,x2) .* problemData.h_tCont(t,x1,x2);

% problemData.fCont{1} = @(t,x1,x2) 0*x1;
% problemData.fCont{2} = @(t,x1,x2) 0*x1;
% problemData.fCont{3} = @(t,x1,x2) 0*x1;

end % function
