% Fills the problemData-struct with all basic configuration options.
% Problem parameters are to be modified inside this routine.

%===============================================================================
%> @file template/configureProblem.m
%>
%> @brief Fills the problemData-struct with all basic configuration options.
%>        Problem parameters are to be modified inside this routine.
%===============================================================================
%>
%> @brief Fills the problemData-struct with all basic configuration options.
%>        Problem parameters are to be modified inside this routine.
%>
%> This routine is called before any other function for the problem.
%> It should define all problem parameters.
%>
%> The struct must provide a positive numeric parameter <code>numSteps</code>
%> that specifies the number of iterative steps (e.g., number of time
%> steps) which are to be executed in the main part of the solution
%> algorithm.
%> For stationary problems set this to 1.
%>
%> @param  problemData  A (probably) empty struct with problem parameters.
%>                      @f$[\text{struct}]@f$
%>
%> @retval problemData  A struct with all necessary parameters and definitions
%>                      for the problem description. @f$[\text{struct}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Hennes Hajduk, Balthasar Reuter, Florian Frank, Vadym Aizinger
%> 
%> @par License
%> @parblock
%> This program is free software: you can redistribute it and/or modify
%> it under the terms of the GNU General Public License as published by
%> the Free Software Foundation, either version 3 of the License, or
%> (at your option) any later version.
%>
%> This program is distributed in the hope that it will be useful,
%> but WITHOUT ANY WARRANTY; without even the implied warranty of
%> MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%> GNU General Public License for more details.
%>
%> You should have received a copy of the GNU General Public License
%> along with this program.  If not, see <http://www.gnu.org/licenses/>.
%> @endparblock
%
function problemData = configureProblem(problemData)

problemData.isHotstart = false;

%% Configuration to use: 
% - 'rotation' calls configureRotation()
% - 'analytical' calls configureAnalyticalTest()
% - 'biological' calls configureBiological()
% - 'ADCIRC' reads 'swe/fort_<name>.15'
problemData = setdefault(problemData, 'configSource', 'ADCIRC');

%% What kind of grid to use:
% - 'square' creates a unit square [0,1]x[0,1] with given pd.hmax,
%   open sea boundary in the east (type 4), and land boundary (type 1) on 
%   all other edges 
% - 'hierarchical' creates a unit square [0,1]x[0,1] with specified hmax
%   and performs uniform refinement according to parameter 'refinement'.
%   Boundary type 4 on east-boundary, 1 on all others.
% - 'ADCIRC' reads grid information from 'swe/fort_<name>.{14,17}'.
problemData = setdefault(problemData, 'gridSource', 'ADCIRC');
problemData = setdefault(problemData, 'refinement', 0);

%% Parameters. 
% Set default values if they are not yet available in problemData
problemData = setdefault(problemData, 'p'         , 1);  % local polynomial degree
problemData = setdefault(problemData, 'ordRK'     , min(problemData.p+1,3));  % order of Runge Kutta time stepper
problemData = setdefault(problemData, 'isVisGrid' , false);  % visualization of grid
problemData = setdefault(problemData, 'maskTol'   , 1.0e-14);  % maximal tolerance of slope for which species are considered constant
problemData = setdefault(problemData, 'isCoupling', false); % Receive velocity coefficients and fluxes from a different model, e.g. 'swe'

%% Parameter check.
assert(problemData.p >= 0 && problemData.p <= 4, 'Polynomial order must be zero to four.')
assert(problemData.ordRK >= 1 && problemData.ordRK <= 3, 'Order of Runge Kutta must be zero to three.')

switch problemData.configSource
  case 'rotation'
    problemData = setdefault(problemData, 'numSpecies', 1);  % number of transported species
  case 'analytical'
    problemData = setdefault(problemData, 'numSpecies', 3);  % number of transported species
  case 'biological'
    problemData = setdefault(problemData, 'numSpecies', 3);  % number of transported species
  case 'ADCIRC'
    problemData = setdefault(problemData, 'numSpecies', 3);  % number of transported species
  otherwise
    error('Invalid config source.')
end % switch

problemData = setdefault(problemData, 'isMask', true(problemData.numSpecies, 1));  % computation only where species is not constant
problemData.maskType = 'vertex-based';

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
    problemData = setdefault(problemData, 'hmax'      , 200);  % maximum edge length of triangle
    problemData = setdefault(problemData, 'numSteps'  , 200);  % number of time steps
    problemData = setdefault(problemData, 'tEnd'      , 500);  % end time
    problemData = configureAnalyticalTest(problemData);
  case 'biological'
    problemData.isSolutionAvailable = false;
    problemData = setdefault(problemData, 'hmax'      , 2^-5);  % maximum edge length of triangle
    problemData = setdefault(problemData, 'numSteps'  , 750);  % number of time steps
    problemData = setdefault(problemData, 'tEnd'      , (problemData.numSteps/3142)*2*pi);  % end time
    problemData = configureBiological(problemData);
  case 'ADCIRC'
    problemData.isSolutionAvailable = false;
    problemData = setdefault(problemData, 'hmax'      , 1);  % maximum edge length of triangle
    problemData = setdefault(problemData, 'numSteps'  , 172800);  % number of time steps
    problemData = setdefault(problemData, 'tEnd'      , 864000);  % end time
    problemData.isSolutionAvailable = false;
    problemData = configureADCIRC(problemData);
  otherwise
    error('Invalid config source.')
end % switch

assert(problemData.hmax > 0, 'Maximum edge length must be positive.')
assert(problemData.numSteps > 0, 'Number of time steps must be positive.')
assert(isequal(~problemData.isMask | cell2mat(problemData.isSlopeLim), ones(problemData.numSpecies,1)), 'Usage of mask relies on slope limiting.');
end % function

%% LeVeque's solid body rotation
function problemData = configureRotation(problemData)

problemData = setdefault(problemData, 'hCont', @(t,x1,x2) 0.1*(x1==x1));
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
  problemData.isSlopeLim{species}  = true; % slope limiting
  problemData.typeSlopeLim{species} = 'linear'; % Type of slope limiter (linear, hierarch_vert, strict)
  
  problemData.outputFrequency{species} = 100; % no visualization of every timestep
  problemData.outputBasename{species}  = ['output' filesep 'solution_' num2str(species) '_' problemData.typeSlopeLim{species}]; % Basename of output files
  problemData.outputTypes{species}     = cellstr('vtk'); % solution output file types
  
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

ds = 0.01; % domainScale
slope = 0.005;
depth = 2;

A = 0.1;
B = 0.1;
C = 0.01;

problemData.zb_xCont = @(x1,x2) slope*(x1==x1);
problemData.zb_yCont = @(x1,x2) slope*(x1==x1);

problemData.hCont = @(t,x1,x2) C*(sin(ds*(x1-t)) + sin(ds*(x2-t))) - slope*(x1+x2) + depth;
problemData.h_tCont = @(t,x1,x2) -ds*C*(cos(ds*(x1-t)) + cos(ds*(x2-t)));
problemData.h_xCont = @(t,x1,x2) ds*C*cos(ds*(x1-t)) - problemData.zb_xCont(x1,x2);
problemData.h_yCont = @(t,x1,x2) ds*C*cos(ds*(x2-t)) - problemData.zb_yCont(x1,x2);

problemData.uCont = @(t,x1,x2) A*sin(ds*(x1-t));
problemData.u_xCont = @(t,x1,x2) ds*A*cos(ds*(x1-t));
problemData.vCont = @(t,x1,x2) B*sin(ds*(x2-t));
problemData.v_yCont = @(t,x1,x2) ds*B*cos(ds*(x2-t));

problemData = setdefault(problemData, 'uHCont', @(t,x1,x2) problemData.hCont(t,x1,x2) .* problemData.uCont(t,x1,x2));
problemData = setdefault(problemData, 'vHCont', @(t,x1,x2) problemData.hCont(t,x1,x2) .* problemData.vCont(t,x1,x2));

% analytical solution
problemData.solCont = {@(t,x1,x2) sin(ds*(x1+x2-t))+2, @(t,x1,x2) cos(ds*(x1+x2-t))+2, @(t,x1,x2) -cos(ds*(x1+x2-t))+2};
sol_tCont = {@(t,x1,x2) -ds*cos(ds*(x1+x2-t)), @(t,x1,x2) ds*sin(ds*(x1+x2-t)), @(t,x1,x2) -ds*sin(ds*(x1+x2-t))};
sol_xCont = {@(t,x1,x2) ds*cos(ds*(x1+x2-t)), @(t,x1,x2) -ds*sin(ds*(x1+x2-t)), @(t,x1,x2) ds*sin(ds*(x1+x2-t))};
sol_yCont = {@(t,x1,x2) ds*cos(ds*(x1+x2-t)), @(t,x1,x2) -ds*sin(ds*(x1+x2-t)), @(t,x1,x2) ds*sin(ds*(x1+x2-t))};

for species = 1:problemData.numSpecies
  problemData.isVisSol{species}    = true; % visualization of solution
  problemData.isSlopeLim{species}  = true; % slope limiting
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
Vm = 0.001;
ks = 1;
Rm = 0.0001;
ep = 0.0001;
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

end % function

%% Biological application
%% LeVeque's solid body rotation
function problemData = configureBiological(problemData)

problemData = setdefault(problemData, 'hCont', @(t,x1,x2) 0.1*(x1==x1));
problemData = setdefault(problemData, 'uCont', @(t,x1,x2) 0.5 - x2);
problemData = setdefault(problemData, 'vCont', @(t,x1,x2) x1 - 0.5);
problemData = setdefault(problemData, 'uHCont', @(t,x1,x2) problemData.hCont(t,x1,x2) .* problemData.uCont(t,x1,x2));
problemData = setdefault(problemData, 'vHCont', @(t,x1,x2) problemData.hCont(t,x1,x2) .* problemData.vCont(t,x1,x2));

%% Coefficients and boundary data (LeVeque's solid body rotation).
G = @(x1, x2, x1_0, x2_0) (1/0.15) * sqrt((x1-x1_0).^2 + (x2-x2_0).^2);
problemData.cH0Cont = { @(x1, x2) ((x1 - 0.5).^2 + (x2 - 0.75).^2 <= 0.0225 & (x1 <= 0.475 | x1 >= 0.525 | x2 >= 0.85)) + ...
                                  (1-G(x1, x2, 0.5, 0.25)) .* ((x1 - 0.5).^2 + (x2 - 0.25).^2 <= 0.0225) + ...
                                  0.25*(1+cos(pi*G(x1, x2, 0.25, 0.5))).*((x1 - 0.25).^2 + (x2 - 0.5).^2 <= 0.0225);
                        @(x1, x2) 2*((x1 - 0.5).^2 + (x2 - 0.75).^2 <= 0.0225 & (x1 > 0.475 & x1 < 0.525 & x2 < 0.85)) + ...
                                  ((x1 - 0.5).^2 + (x2 - 0.25).^2 <= 0.035) .* ((x1 - 0.5).^2 + (x2 - 0.25).^2 > 0.0225) + ...
                                  ((x1 - 0.25).^2 + (x2 - 0.5).^2 <= 0.02); 
                        @(x1, x2) -cos(2*pi*(x1-0.5)) .* sin(2*pi*x2) .* (x1 > 0.25 & x1 < 0.75 & x2 > 0.5)};

for species = 1:problemData.numSpecies
  problemData.isVisSol{species}    = true; % visualization of solution
  problemData.isSlopeLim{species}  = true; % slope limiting
  problemData.typeSlopeLim{species} = 'linear'; % Type of slope limiter (linear, hierarch_vert, strict)
  
  problemData.outputFrequency{species} = 100; % no visualization of every timestep
  problemData.outputBasename{species}  = ['output' filesep 'solution_' num2str(species) '_' problemData.typeSlopeLim{species}]; % Basename of output files
  problemData.outputTypes{species}     = cellstr('vtk'); % solution output file types
  
  %% Parameter check.
  assert(~problemData.isSlopeLim{species} || problemData.p > 0, 'Slope limiting only available for p > 0.')

  problemData.cDCont{species} = @(t,x1,x2) zeros(size(x1));
  problemData.gNCont{species} = @(t,x1,x2) zeros(size(x1));
end % for

%% Reaction term definitions.
% NPZ model
% parameters (Notation as in paper)
I0 = 1;
Vm = 0.001;
ks = 1;
Rm = 0.1;
ep = 0.001;
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

problemData.fCont{1} = @(t,x1,x2) 0*x1;
problemData.fCont{2} = @(t,x1,x2) 0*x1;
problemData.fCont{3} = @(t,x1,x2) 0*x1;

end % function

%% Coupling to ADCIRC grid
function problemData = configureADCIRC(problemData)

problemData = setdefault(problemData, 'name', 'galv');
numForcings = 0; % since this information is unnecessary at this point, it will not be computed
isSpherical = 0;
projCenter = [0 0];

problemData = setdefault(problemData, 'hCont', @(t,x1,x2) x1==x1);
problemData = setdefault(problemData, 'uHCont', @(t,x1,x2) 0*x1);
problemData = setdefault(problemData, 'vHCont', @(t,x1,x2) 0*x1);

problemData = setdefault(problemData, 'domainADCIRC', getFunctionHandle('swe/domainADCIRC'));

[problemData.g, depth] = problemData.domainADCIRC(['swe/fort_' problemData.name '.14'], ['swe/fort_' problemData.name '.17'], numForcings, isSpherical, projCenter);

problemData.isHotstart = true;
hotstartData = readHotstart('output/galv_1.mat');

zbCont = @(x1,x2) execin('swe/evaluateFuncFromVertexValues', problemData.g, -depth, x1, x2);

N = nchoosek(problemData.p+2, problemData.p);
% TODO because of zb this will not be consistent to swe
problemData.h0Disc = hotstartData.cDisc(:,:,1) - projectFuncCont2DataDisc(problemData.g, zbCont, 2*problemData.p+1, eye(N), computeBasesOnQuad(N, struct));

problemData.xiOSCont = @(t,x1,x2) ( cos(0.000067597751162*t) * 0.075 * cos(-194.806 * pi/180) ...
                                  - sin(0.000067597751162*t) * 0.075 * sin(-194.806 * pi/180) ...
                                  + cos(0.000072921165921*t) * 0.095 * cos(-206.265 * pi/180) ...
                                  - sin(0.000072921165921*t) * 0.095 * sin(-206.265 * pi/180) ...
                                  + cos(0.000137879713787*t) * 0.100 * cos(-340.000 * pi/180) ...
                                  - sin(0.000137879713787*t) * 0.100 * sin(-340.000 * pi/180) ...
                                  + cos(0.000140518917083*t) * 0.395  ... % ... * cos(0) - ... * sin(0)
                                  + cos(0.000145444119418*t) * 0.060 * cos(-42.9718 * pi/180) ...
                                  - sin(0.000145444119418*t) * 0.060 * sin(-42.9718 * pi/180) ) * (x2 <= 3280000) + (x2 > 3280000);

aux = false(3397,N,2);
aux([725 726 727 794 795 796 797 870 871 872 933 934 997 998 999 1000 1585 1586 1587 1625 1626 1627 1628 1629 1670 1671 1672 1673],:,1) = true;
aux([423 424 465 466 467 468 805 806 807 876 877 878 1366 1367 1368 1424 1425 1426 1427],:,2) = true;
problemData.cH0Disc = { 1E-6 * aux(:,:,1) .* problemData.h0Disc;
                        1E-6 * aux(:,:,2) .* problemData.h0Disc;
                        zeros(3397,N) };

problemData.cDCont = { @(t,x1,x2) 0*x1; @(t,x1,x2) 0*x1; @(t,x1,x2) 7.0*1E-6*(x2 > 3280000) }; % besprechen, 14 g/mol * 0.5 mmol/m^3 = 7.0*1E-6 kg / m^3

problemData.outputBasename = {['output' filesep 'phyto'], ['output' filesep 'zoo'], ['output' filesep 'nitro']};

for species = 1:problemData.numSpecies
  problemData.isVisSol{species}    = true; % visualization of solution
  problemData.isSlopeLim{species}  = true; % slope limiting
  problemData.typeSlopeLim{species} = 'linear'; % Type of slope limiter (linear, hierarch_vert, strict)
  
  problemData.outputFrequency{species} = 8640; % no visualization of every timestep
  problemData.outputTypes{species}     = cellstr('vtk'); % solution output file types
  
  %% Parameter check.
  assert(~problemData.isSlopeLim{species} || problemData.p > 0, 'Slope limiting only available for p > 0.')

  problemData.gNCont{species} = @(t,x1,x2) zeros(size(x1));
end % for

%% Reaction term definitions.
% NPZ model
% parameters (Notation as in paper)
I0 = 35; % muE / (m s^2)
Vm = 2.0 / 86400; % 1 / s
ks = 2.0E-7; % 0.2 mug / l = 0.2E-9 / (0.1m)^3 = 2.0E-7 kg / m^3
Rm = 0.5 / 86400; % 1 / s
ep1 = 0.1 / 86400; % 1 / s
ep2 = 0.2 / 86400; % 1 / s
gamma = 0.7; % no unit
lambda = 5.0E5; % 0.5 (0.1m)^3 / (mug) = 5.0 1E-4 / 1E-9 m^3 / kg
h0 = 10; % m
k = 0.1; % m

I = @(t,x1,x2) I0 * h0 / k * (1 - exp(k*zbCont(x1,x2))) / (1 - exp(-k*h0));
f = @(t,x1,x2) 1 - exp(I(t,x1,x2) / I0); % saturating response
g = @(t,x1,x2,N) Vm * N ./ (ks + N); % Michaelis-Menten uptake
h = @(t,x1,x2,P) Rm * (1 - exp(-lambda * P)); % saturating (Ivlev)
i = @(t,x1,x2,P) ep1; % linear death rate
j = @(t,x1,x2,Z) ep2; % linear death rate

problemData.reactions{1} = @(t,x1,x2,c,cH) f(t,x1,x2) .* g(t,x1,x2,c{3}) .* cH{1} - h(t,x1,x2,c{1}) .* cH{2} - i(t,x1,x2,c{1}) .* cH{1};
problemData.reactions{2} = @(t,x1,x2,c,cH) gamma * h(t,x1,x2,c{1}) .* cH{2} - j(t,x1,x2,c{2}) .* cH{2};
problemData.reactions{3} = @(t,x1,x2,c,cH) -f(t,x1,x2) .* g(t,x1,x2,c{3}) .* cH{1} + (1-gamma) .* h(t,x1,x2,c{1}) .* cH{2} + i(t,x1,x2,c{1}) .* cH{1} + j(t,x1,x2,c{2}) .* cH{2};

problemData.fCont{1} = @(t,x1,x2) 0*x1;
problemData.fCont{2} = @(t,x1,x2) 0*x1;
problemData.fCont{3} = @(t,x1,x2) 0*x1;

end % function
