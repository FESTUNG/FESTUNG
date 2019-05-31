% Fills the problemData-struct with all basic configuration options.
% Problem parameters are to be modified inside this routine.

%===============================================================================
%> @file ./cahn-hilliard/configureProblem.m
%>
%> @brief Fills the problemData-struct with all basic configuration options.
%>        Problem parameters are to be modified inside this routine.
%===============================================================================
%>
%> @brief Fills the problemData-struct with all basic configuration options.
%>        Problem parameters are to be modified inside this routine.
%>
%> This routine is called before any other function for the problem.
%> It defines all problem parameters and should be the only file users have
%> to adapt to their needs.
%>
%> It provides all configuration options for the numerical solver to
%> approximate solutions
%> @f$c:\overline{J}\times\overline{\Omega}\rightarrow\mathbb{R}@f$ and
%> @f$\Phi: \overline{J}\times\overline{\Omega}\rightarrow\mathbb{R}@f$
%> of the Cahn-Hilliard equation
%> @f{align*}{
%> \partial_t c + \nabla \cdot m \nabla \Phi &\;=\; f                 &&\text{in}~J\times\Omega\,,\\
%> \Psi'(c) - \epsilon^2 \Delta C            &\;=\; \Phi              &&\text{in}~J\times\Omega\,,\\
%> \nabla c \cdot \nu                        &\;=\; g_\mathrm{N, con} &&\text{on}~J\times\partial\Omega\,,\\
%> \nabla \Phi \cdot \nu                     &\;=\; g_\mathrm{N, pot} &&\text{on}~J\times\partial\Omega\,,\\
%> c                                         &\;=\; c^0               &&\text{on}~\{0\}\times\Omega\,.
%> @f}
%> The right hand side@f$f:J\times\Omega\rightarrow \mathbb{R}@f$ and boundary
%> data @f$g_\mathrm{N, con}, \, g_\mathrm{N, pot}@f$ may vary in time and space.
%> A detailed description can be found in @ref RAWFK2016.
%>
%> Please read the inline-comments in the code for the meaning of each
%> configuration option.
%>
%> @param  problemData  A (probably) empty struct with problem parameters.
%>                      @f$[\text{struct}]@f$
%>
%> @retval problemData  A struct with all necessary parameters and definitions
%>                      for the problem description. @f$[\text{struct}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Balthasar Reuter, Florian Frank, Vadym Aizinger
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
%% Parameters.
% Maximum edge length of triangle
problemData = setdefault(problemData, 'hmax', 2^-6);

% Define simulated scenario ('conv-test', 'rand', 'ball', 'merge-balls')
problemData = setdefault(problemData, 'scenario', 'rand');
% Additional parameter for conv-test-scenario('const', 'cos', 'deg', 'deg2', 'exp')
problemData = setdefault(problemData,'conv_test_mobility', 'const');

% Define which potential should be used ('Ginzburg-Landau', 'logarithmic')
problemData = setdefault(problemData, 'potential', 'double_obstacle');

% Configure the DG method
problemData = setdefault(problemData, 'p', 1); % Local polynomial approximation order (0 to 4)
problemData = setdefault(problemData, 'eta', 1); % -1 = SIPG, 0 = IIPG, 1 = NIPG
problemData = setdefault(problemData, 'sigma', 1); % Penalty parameter

% Time stepping parameters
problemData = setdefault(problemData, 'tau', 1e-4); % time step size
problemData = setdefault(problemData, 'tEnd', 1); % end time
problemData = setdefault(problemData, 'implicit', true); % implicit or explicit time stepping

% Adaptive time stepping
problemData = setdefault(problemData, 'adaptiveStepping', false); % activate adaptive time stepping
problemData = setdefault(problemData, 'tauMax', 1e-1); % maximum time step size
problemData = setdefault(problemData, 'tauMin', 1e-8); % minimum time step size
problemData = setdefault(problemData, 'tauDecr', 0.25); % factor time step size is reduced
problemData = setdefault(problemData, 'tauIncr', 2); % factor time step size in enlarged
problemData = setdefault(problemData, 'newtonStepRatio', 0.25); % if less than newtonStepRatio * maxNewtonSteps steps used enlarge

% Slope limiting settings
problemData = setdefault(problemData, 'isFluxLim', false); % enable/disable slope limiting
problemData = setdefault(problemData, 'gamma', 5e0); % Limiting parameter
problemData = setdefault(problemData, 'beta', 1e-4); % Regularization parameter beta > 0 (small)
problemData = setdefault(problemData, 'initialLimiting', true); % enable/disable slope limiting of initial data
problemData = setdefault(problemData, 'typeSlopeLim', 'hierarch_vert'); % type of initial limitation
problemData = setdefault(problemData, 'standardLim', true); % Decide whether standard slope limiting technique should be used
problemData = setdefault(problemData, 'zalesak', true); % Decide whether Zalesak's algorithm should be applied as post-processing
problemData = setdefault(problemData, 'zalesakSteps', 100); % Number of Iterations within Zalesak's postprocessing scheme

% Visualization settings
problemData = setdefault(problemData, 'isVisGrid', false); % visualization of grid
problemData = setdefault(problemData, 'isVisSol', true); % visualization of solution
problemData = setdefault(problemData, 'outputFrequency', 10); % no visualization of every timestep
problemData = setdefault(problemData, 'outputBasename', ...
                         ['output' filesep 'cahnHilliard']); % Define name of output file
problemData = setdefault(problemData, 'outputTypes', { 'vtk' }); % Type of visualization files ('vtk, 'tec')

% Decide whether Voronoi or standard triangulation should be used
problemData = setdefault(problemData, 'Voronoi_tri', false); % default should be false
% Decide whether Friedrich-Keller triangulation should be used
problemData = setdefault(problemData, 'Friedrichs_Keller', false);

% Decide whether MATLAB should inform about the simulation being finished
problemData = setdefault(problemData, 'gong', false);

% Newton-Krylow Method Configuration(1,2,3)
% option 1: The mobility-dependent matrix is evaluated at C_old the solution of the old time-step
% option 2: The mobility-dependent matrix is updated after each newton-step
% option 3: The mobility-dependent matrix is treated implicitly with respect to time. The non-linear equation is solved via newton-krylov
% option 4: The mobility-dependent matrix is treated balancedly between opt. 3 and opt 4 (default)
problemData = setdefault(problemData, 'newton_version', 4);
% Preconditioner for newton-krylov,newton_version = 3 (1,2)
% option 1: no preconditioner (is currently implemented via P = eye)
% option 2: uses jacobian of linear part (where no mobility dependent parts has been evaluated)
% option 3: uses jacobian of full problem
problemData = setdefault(problemData, 'newton_precond', 3);

%% Coefficients and boundary data
if strcmp(problemData.potential,'Ginzburg-Landau')
  problemData.ConvexPsi  = @(C) C.^3; % Derivative of the Convex Part of Potential Psi
  problemData.ConcavePsi = @(C) -C; % Derivative of the Concave Part of Potential Psi
  problemData.DerivativeOfPsi = @(C) problemData.ConvexPsi(C) + problemData.ConcavePsi(C);
  problemData.JacobiConvexPsi = @(C) 3 * C.^2; % Jacobian of the implicit part of Psi
elseif strcmp(problemData.potential,'logarithmic')
  vartheta = 0.4; % absolute temperature
  vartheta_c = 1; % critical temperature (bigger than temperature for double well potential)
  problemData.ConvexPsi  = @(C) vartheta * 0.5 * ( log(1+C) - log(1-C) ) ; % Derivative of the Convex Part of Potential Psi
  problemData.ConcavePsi = @(C) - vartheta_c * C; % Derivative of the Concave Part of Potential Psi
  problemData.JacobiConvexPsi = @(C) - vartheta ./ (C.^2 - 1); % Jacobian of the implicit part of Psi
elseif strcmp(problemData.potential,'double_obstacle')
  problemData.ConvexPsi  = @(C) zeros(size(C)); % Derivative of the Convex Part of Potential Psi
  problemData.ConcavePsi = @(C) obstaclePotentialPsiPrime(C);
  problemData.DerivativeOfPsi = @(C) problemData.ConvexPsi(C) + problemData.ConcavePsi(C);
  problemData.JacobiConvexPsi = @(C) zeros(size(C)); % Jacobian of the implicit part of Psi
end
problemData.mobilityCont =  @(Y) ones(size(Y));
if strcmp(problemData.scenario,'conv-test')
  % Convergence test for Cahn-Hilliard
  switch problemData.conv_test_mobility
    case 'const'
      problemData.epsilon = 1;
      problemData.mobility = 1;
      if problemData.Voronoi_tri
        problemData.cCont = @(t,x,y) exp(-t) * sin(2 * pi * (x - y / sqrt(3))) .* sin(4 * pi * y / sqrt(3));
        problemData.c0Cont = @(x,y) problemData.cCont(0,x,y);
        problemData.fCont = @(t,x,y) (exp(-3*t)*(9*exp(2*t)*cos((2*pi*(3*x + 3^(1/2)*y))/3) ...
          + 81*pi^2*cos(2*pi*(x - 3^(1/2)*y)) + 81*pi^2*cos(6*pi*(x - 3^(1/2)*y)) ...
          - 27*pi^2*cos(2*pi*(3*x + 3^(1/2)*y)) - 27*pi^2*cos((2*pi*(3*x + 3^(1/2)*y))/3) ...
          + 63*pi^2*cos((2*pi*(3*x + 5*3^(1/2)*y))/3) - 117*pi^2*cos((2*pi*(3*x - 7*3^(1/2)*y))/3) ...
          + 63*pi^2*cos((2*pi*(9*x - 3^(1/2)*y))/3) - 117*pi^2*cos((2*pi*(9*x - 5*3^(1/2)*y))/3) ...
          - 9*cos(2*pi*(x - 3^(1/2)*y))*exp(2*t) - 144*pi^2*cos(2*pi*(x - 3^(1/2)*y))*exp(2*t) ...
          + 2304*pi^4*cos(2*pi*(x - 3^(1/2)*y))*exp(2*t) + 48*pi^2*exp(2*t)*cos((2*pi*(3*x + 3^(1/2)*y))/3) ...
          - 256*pi^4*exp(2*t)*cos((2*pi*(3*x + 3^(1/2)*y))/3)))/18;
        grad_c = @(t,x,y) [ 2*pi*cos(2*pi*(x - (3^(1/2)*y)/3))*exp(-t).*sin((4*pi*3^(1/2)*y)/3) , ...
          (4*3^(1/2)*pi*sin(2*pi*(x - (3^(1/2)*y)/3))*exp(-t).*cos((4*pi*3^(1/2)*y)/3))/3 ...
          - (2*3^(1/2)*pi*cos(2*pi*(x - (3^(1/2)*y)/3))*exp(-t).*sin((4*pi*3^(1/2)*y)/3))/3 ];
        problemData.gNContCon = @(t,x,y) ...
          + grad_c(t,x,y) * [0; 1; 0; 1] .* (y == sqrt(3)/2) + grad_c(t,x,y) * [0; -1; 0; -1] .* (y == 0) ...
          + grad_c(t,x,y) * [-sqrt(3); 1; -sqrt(3); 1] ./ 2 .* ((x < 0.8) .* (y > 0) .* (y < sqrt(3)/2)) ...
          + grad_c(t,x,y) * [sqrt(3); -1; sqrt(3); -1] ./ 2 .* ((x > 0.8) .* (y > 0) .* (y < sqrt(3)/2));
        grad_mu = @(t,x,y) [(64*pi^3*cos(2*pi*(x - (3^(1/2)*y)/3))*exp(-t).*sin((4*pi*3^(1/2)*y)/3))/3 ...
          - (32*pi^3*sin(2*pi*(x - (3^(1/2)*y)/3))*exp(-t).*cos((4*pi*3^(1/2)*y)/3))/3 ...
          + 2*pi*cos(2*pi*(x - (3^(1/2)*y)/3))*exp(-t).*sin((4*pi*3^(1/2)*y)/3).*(sin(2*pi*(x - (3^(1/2)*y)/3)).^2 ...
          *exp(-2*t).*sin((4*pi*3^(1/2)*y)/3).^2 - 1) ...
          + 4*pi*cos(2*pi*(x - (3^(1/2)*y)/3)).*sin(2*pi*(x - (3^(1/2)*y)/3)).^2*exp(-3*t).*sin((4*pi*3^(1/2)*y)/3).^3 , ...
          sin(2*pi*(x - (3^(1/2)*y)/3))*exp(-t).*sin((4*pi*3^(1/2)*y)/3).*((8*3^(1/2)*pi*sin(2*pi*(x - (3^(1/2)*y)/3)).^2 ...
          *exp(-2*t).*cos((4*pi*3^(1/2)*y)/3).*sin((4*pi*3^(1/2)*y)/3))/3 - (4*3^(1/2)*pi*cos(2*pi*(x - (3^(1/2)*y)/3)) ...
          .*sin(2*pi*(x - (3^(1/2)*y)/3))*exp(-2*t).*sin((4*pi*3^(1/2)*y)/3).^2)/3) - (128*3^(1/2)*pi^3*cos(2*pi*(x - (3^(1/2)*y)/3)) ...
          *exp(-t).*sin((4*pi*3^(1/2)*y)/3))/9 + (160*3^(1/2)*pi^3*sin(2*pi*(x - (3^(1/2)*y)/3))*exp(-t).*cos((4*pi*3^(1/2)*y)/3))/9  ...
          - (2*3^(1/2)*pi*cos(2*pi*(x - (3^(1/2)*y)/3))*exp(-t).*sin((4*pi*3^(1/2)*y)/3).*(sin(2*pi*(x - (3^(1/2)*y)/3)).^2*exp(-2*t) ...
          .*sin((4*pi*3^(1/2)*y)/3).^2 - 1))/3 + (4*3^(1/2)*pi*sin(2*pi*(x - (3^(1/2)*y)/3))*exp(-t)...
          .*cos((4*pi*3^(1/2)*y)/3).*(sin(2*pi*(x - (3^(1/2)*y)/3)).^2*exp(-2*t).*sin((4*pi*3^(1/2)*y)/3).^2 - 1))/3];
        problemData.gNContPot = @(t,x,y) ...
          + grad_mu(t,x,y) * [0; 1; 0 ;1] .* (y == sqrt(3)/2) + grad_mu(t,x,y) * [0; -1; 0; -1] .* (y == 0) ...
          + grad_mu(t,x,y) * [-sqrt(3); 1; -sqrt(3); 1] ./ 2 .* ((x < 0.8) .* (y > 0) .* (y < sqrt(3)/2)) ...
          + grad_mu(t,x,y) * [sqrt(3); -1; sqrt(3); -1] ./ 2 .* ((x > 0.8) .* (y > 0) .* (y < sqrt(3)/2));
      else
        problemData.cCont = @(t,x,y) exp(-t) * sin(2 * pi * x) .* sin(2 * pi * y);
        problemData.c0Cont = @(x,y) problemData.cCont(0,x,y);
        problemData.fCont = @(t,x,y) -exp(-3*t) * sin(2*pi*x) .* sin(2*pi*y) .* ...
          (exp(2*t) + 24 * pi^2 * sin(2*pi*x).^2 + 24 * pi^2 * sin(2*pi*y).^2 + 8 * pi^2 * exp(2*t) ...
          - 64 * pi^4 * exp(2*t) - 72 * pi^2 * sin(2*pi*x).^2 .* sin(2*pi*y).^2);
        problemData.gNContCon = @(t,x,y) ...
          - 2 * pi * exp(-t) * sin(2*pi*y) .* (x == 0) ...
          + 2 * pi * exp(-t) * sin(2*pi*y) .* (x == 1) ...
          - 2 * pi * exp(-t) * sin(2*pi*x) .* (y == 0) ...
          + 2 * pi * exp(-t) * sin(2*pi*x) .* (y == 1);
        problemData.gNContPot = @(t,x,y) ...
          + (2*pi*exp(-t)*sin(2*pi*y) - 16*pi^3*exp(-t)*sin(2*pi*y)) .* (x == 0) ...
          + (16*pi^3*exp(-t)*sin(2*pi*y) - 2*pi*exp(-t)*sin(2*pi*y)) .* (x == 1) ...
          + (2*pi*exp(-t)*sin(2*pi*x) - 16*pi^3*exp(-t)*sin(2*pi*x)) .* (y == 0) ...
          + (16*pi^3*exp(-t)*sin(2*pi*x) - 2*pi*exp(-t)*sin(2*pi*x)) .* (y == 1);
      end
    case 'cos'
      problemData.epsilon = 1;
      problemData.mobility = 'non-constant';
      problemData.mobilityCont = @(c) cos(c);
      problemData.cCont = @(t,x,y) exp(-t) * sin(2 * pi * x) .* sin(2 * pi * y);
      problemData.c0Cont = @(x,y) problemData.cCont(0,x,y);
      problemData.fCont = @(t,x,y)cos(exp(-t).*sin(x.*pi.*2.0).*sin(y.*pi.*2.0)).*(pi.^2.*exp(t.*-3.0).*sin(x.*pi.*2.0).^3.*sin(y.*pi.*2.0).^3.*8.0+pi.^4.*exp(-t).*sin(x.*pi.*2.0).*sin(y.*pi.*2.0).*3.2e1+pi.^2.*exp(-t).*sin(x.*pi.*2.0).*sin(y.*pi.*2.0).*(exp(t.*-2.0).*sin(x.*pi.*2.0).^2.*sin(y.*pi.*2.0).^2-1.0).*4.0-pi.^2.*exp(t.*-3.0).*cos(x.*pi.*2.0).^2.*sin(x.*pi.*2.0).*sin(y.*pi.*2.0).^3.*2.4e1)+cos(exp(-t).*sin(x.*pi.*2.0).*sin(y.*pi.*2.0)).*(pi.^2.*exp(t.*-3.0).*sin(x.*pi.*2.0).^3.*sin(y.*pi.*2.0).^3.*8.0+pi.^4.*exp(-t).*sin(x.*pi.*2.0).*sin(y.*pi.*2.0).*3.2e1+pi.^2.*exp(-t).*sin(x.*pi.*2.0).*sin(y.*pi.*2.0).*(exp(t.*-2.0).*sin(x.*pi.*2.0).^2.*sin(y.*pi.*2.0).^2-1.0).*4.0-pi.^2.*exp(t.*-3.0).*cos(y.*pi.*2.0).^2.*sin(x.*pi.*2.0).^3.*sin(y.*pi.*2.0).*2.4e1)-exp(-t).*sin(x.*pi.*2.0).*sin(y.*pi.*2.0)+pi.*sin(exp(-t).*sin(x.*pi.*2.0).*sin(y.*pi.*2.0)).*exp(-t).*cos(x.*pi.*2.0).*sin(y.*pi.*2.0).*(pi.^3.*exp(-t).*cos(x.*pi.*2.0).*sin(y.*pi.*2.0).*1.6e1+pi.*exp(t.*-3.0).*cos(x.*pi.*2.0).*sin(x.*pi.*2.0).^2.*sin(y.*pi.*2.0).^3.*4.0+pi.*exp(-t).*cos(x.*pi.*2.0).*sin(y.*pi.*2.0).*(exp(t.*-2.0).*sin(x.*pi.*2.0).^2.*sin(y.*pi.*2.0).^2-1.0).*2.0).*2.0+pi.*sin(exp(-t).*sin(x.*pi.*2.0).*sin(y.*pi.*2.0)).*exp(-t).*cos(y.*pi.*2.0).*sin(x.*pi.*2.0).*(pi.^3.*exp(-t).*cos(y.*pi.*2.0).*sin(x.*pi.*2.0).*1.6e1+pi.*exp(t.*-3.0).*cos(y.*pi.*2.0).*sin(x.*pi.*2.0).^3.*sin(y.*pi.*2.0).^2.*4.0+pi.*exp(-t).*cos(y.*pi.*2.0).*sin(x.*pi.*2.0).*(exp(t.*-2.0).*sin(x.*pi.*2.0).^2.*sin(y.*pi.*2.0).^2-1.0).*2.0).*2.0;
      problemData.gNContCon = @(t,x,y) ...
        - 2 * pi * exp(-t) * sin(2*pi*y) .* (x == 0) ...
        + 2 * pi * exp(-t) * sin(2*pi*y) .* (x == 1) ...
        - 2 * pi * exp(-t) * sin(2*pi*x) .* (y == 0) ...
        + 2 * pi * exp(-t) * sin(2*pi*x) .* (y == 1);
      problemData.gNContPot = @(t,x,y) ...
        + (2*pi*exp(-t)*sin(2*pi*y) - 16*pi^3*exp(-t)*sin(2*pi*y)) .* (x == 0) ...
        + (16*pi^3*exp(-t)*sin(2*pi*y) - 2*pi*exp(-t)*sin(2*pi*y)) .* (x == 1) ...
        + (2*pi*exp(-t)*sin(2*pi*x) - 16*pi^3*exp(-t)*sin(2*pi*x)) .* (y == 0) ...
        + (16*pi^3*exp(-t)*sin(2*pi*x) - 2*pi*exp(-t)*sin(2*pi*x)) .* (y == 1);
    case 'deg'
      % increase penalty parameter!!!
      problemData.epsilon = 1;
      problemData.mobility = 'non-constant';
      problemData.mobilityCont = @(c) max(1 - c.^2,0) ;
      problemData.cCont = @(t,x,y) exp(-t) * sin(2 * pi * x) .* sin(2 * pi * y);
      problemData.c0Cont = @(x,y) problemData.cCont(0,x,y);
      problemData.fCont = @(t,x,y)-(exp(t.*-2.0).*sin(x.*pi.*2.0).^2.*sin(y.*pi.*2.0).^2-1.0).*(pi.^2.*exp(t.*-3.0).*sin(x.*pi.*2.0).^3.*sin(y.*pi.*2.0).^3.*8.0+pi.^4.*exp(-t).*sin(x.*pi.*2.0).*sin(y.*pi.*2.0).*3.2e1+pi.^2.*exp(-t).*sin(x.*pi.*2.0).*sin(y.*pi.*2.0).*(exp(t.*-2.0).*sin(x.*pi.*2.0).^2.*sin(y.*pi.*2.0).^2-1.0).*4.0-pi.^2.*exp(t.*-3.0).*cos(x.*pi.*2.0).^2.*sin(x.*pi.*2.0).*sin(y.*pi.*2.0).^3.*2.4e1)-(exp(t.*-2.0).*sin(x.*pi.*2.0).^2.*sin(y.*pi.*2.0).^2-1.0).*(pi.^2.*exp(t.*-3.0).*sin(x.*pi.*2.0).^3.*sin(y.*pi.*2.0).^3.*8.0+pi.^4.*exp(-t).*sin(x.*pi.*2.0).*sin(y.*pi.*2.0).*3.2e1+pi.^2.*exp(-t).*sin(x.*pi.*2.0).*sin(y.*pi.*2.0).*(exp(t.*-2.0).*sin(x.*pi.*2.0).^2.*sin(y.*pi.*2.0).^2-1.0).*4.0-pi.^2.*exp(t.*-3.0).*cos(y.*pi.*2.0).^2.*sin(x.*pi.*2.0).^3.*sin(y.*pi.*2.0).*2.4e1)-exp(-t).*sin(x.*pi.*2.0).*sin(y.*pi.*2.0)+pi.*exp(t.*-2.0).*cos(x.*pi.*2.0).*sin(x.*pi.*2.0).*sin(y.*pi.*2.0).^2.*(pi.^3.*exp(-t).*cos(x.*pi.*2.0).*sin(y.*pi.*2.0).*1.6e1+pi.*exp(t.*-3.0).*cos(x.*pi.*2.0).*sin(x.*pi.*2.0).^2.*sin(y.*pi.*2.0).^3.*4.0+pi.*exp(-t).*cos(x.*pi.*2.0).*sin(y.*pi.*2.0).*(exp(t.*-2.0).*sin(x.*pi.*2.0).^2.*sin(y.*pi.*2.0).^2-1.0).*2.0).*4.0+pi.*exp(t.*-2.0).*cos(y.*pi.*2.0).*sin(x.*pi.*2.0).^2.*sin(y.*pi.*2.0).*(pi.^3.*exp(-t).*cos(y.*pi.*2.0).*sin(x.*pi.*2.0).*1.6e1+pi.*exp(t.*-3.0).*cos(y.*pi.*2.0).*sin(x.*pi.*2.0).^3.*sin(y.*pi.*2.0).^2.*4.0+pi.*exp(-t).*cos(y.*pi.*2.0).*sin(x.*pi.*2.0).*(exp(t.*-2.0).*sin(x.*pi.*2.0).^2.*sin(y.*pi.*2.0).^2-1.0).*2.0).*4.0;
      problemData.gNContCon = @(t,x,y) ...
        - 2 * pi * exp(-t) * sin(2*pi*y) .* (x == 0) ...
        + 2 * pi * exp(-t) * sin(2*pi*y) .* (x == 1) ...
        - 2 * pi * exp(-t) * sin(2*pi*x) .* (y == 0) ...
        + 2 * pi * exp(-t) * sin(2*pi*x) .* (y == 1);
      problemData.gNContPot = @(t,x,y) ...
        + (2*pi*exp(-t)*sin(2*pi*y) - 16*pi^3*exp(-t)*sin(2*pi*y)) .* (x == 0) ...
        + (16*pi^3*exp(-t)*sin(2*pi*y) - 2*pi*exp(-t)*sin(2*pi*y)) .* (x == 1) ...
        + (2*pi*exp(-t)*sin(2*pi*x) - 16*pi^3*exp(-t)*sin(2*pi*x)) .* (y == 0) ...
        + (16*pi^3*exp(-t)*sin(2*pi*x) - 2*pi*exp(-t)*sin(2*pi*x)) .* (y == 1);
    case 'deg2'
      % increase penalty parameter!!!
      problemData.epsilon = 1;
      problemData.mobility = 'non-constant';
      problemData.mobilityCont = @(c) 1 - 0.9 * c.^2 ;
      problemData.cCont = @(t,x,y) exp(-t) * sin(2 * pi * x) .* sin(2 * pi * y);
      problemData.c0Cont = @(x,y) problemData.cCont(0,x,y);
      problemData.fCont = @(t,x,y)-(exp(t.*-2.0).*sin(x.*pi.*2.0).^2.*sin(y.*pi.*2.0).^2.*(9.0./1.0e1)-1.0).*(pi.^2.*exp(t.*-3.0).*sin(x.*pi.*2.0).^3.*sin(y.*pi.*2.0).^3.*8.0+pi.^4.*exp(-t).*sin(x.*pi.*2.0).*sin(y.*pi.*2.0).*3.2e1+pi.^2.*exp(-t).*sin(x.*pi.*2.0).*sin(y.*pi.*2.0).*(exp(t.*-2.0).*sin(x.*pi.*2.0).^2.*sin(y.*pi.*2.0).^2-1.0).*4.0-pi.^2.*exp(t.*-3.0).*cos(x.*pi.*2.0).^2.*sin(x.*pi.*2.0).*sin(y.*pi.*2.0).^3.*2.4e1)-(exp(t.*-2.0).*sin(x.*pi.*2.0).^2.*sin(y.*pi.*2.0).^2.*(9.0./1.0e1)-1.0).*(pi.^2.*exp(t.*-3.0).*sin(x.*pi.*2.0).^3.*sin(y.*pi.*2.0).^3.*8.0+pi.^4.*exp(-t).*sin(x.*pi.*2.0).*sin(y.*pi.*2.0).*3.2e1+pi.^2.*exp(-t).*sin(x.*pi.*2.0).*sin(y.*pi.*2.0).*(exp(t.*-2.0).*sin(x.*pi.*2.0).^2.*sin(y.*pi.*2.0).^2-1.0).*4.0-pi.^2.*exp(t.*-3.0).*cos(y.*pi.*2.0).^2.*sin(x.*pi.*2.0).^3.*sin(y.*pi.*2.0).*2.4e1)-exp(-t).*sin(x.*pi.*2.0).*sin(y.*pi.*2.0)+pi.*exp(t.*-2.0).*cos(x.*pi.*2.0).*sin(x.*pi.*2.0).*sin(y.*pi.*2.0).^2.*(pi.^3.*exp(-t).*cos(x.*pi.*2.0).*sin(y.*pi.*2.0).*1.6e1+pi.*exp(t.*-3.0).*cos(x.*pi.*2.0).*sin(x.*pi.*2.0).^2.*sin(y.*pi.*2.0).^3.*4.0+pi.*exp(-t).*cos(x.*pi.*2.0).*sin(y.*pi.*2.0).*(exp(t.*-2.0).*sin(x.*pi.*2.0).^2.*sin(y.*pi.*2.0).^2-1.0).*2.0).*(1.8e1./5.0)+pi.*exp(t.*-2.0).*cos(y.*pi.*2.0).*sin(x.*pi.*2.0).^2.*sin(y.*pi.*2.0).*(pi.^3.*exp(-t).*cos(y.*pi.*2.0).*sin(x.*pi.*2.0).*1.6e1+pi.*exp(t.*-3.0).*cos(y.*pi.*2.0).*sin(x.*pi.*2.0).^3.*sin(y.*pi.*2.0).^2.*4.0+pi.*exp(-t).*cos(y.*pi.*2.0).*sin(x.*pi.*2.0).*(exp(t.*-2.0).*sin(x.*pi.*2.0).^2.*sin(y.*pi.*2.0).^2-1.0).*2.0).*(1.8e1./5.0);
      problemData.gNContCon = @(t,x,y) ...
        - 2 * pi * exp(-t) * sin(2*pi*y) .* (x == 0) ...
        + 2 * pi * exp(-t) * sin(2*pi*y) .* (x == 1) ...
        - 2 * pi * exp(-t) * sin(2*pi*x) .* (y == 0) ...
        + 2 * pi * exp(-t) * sin(2*pi*x) .* (y == 1);
      problemData.gNContPot = @(t,x,y) ...
        + (2*pi*exp(-t)*sin(2*pi*y) - 16*pi^3*exp(-t)*sin(2*pi*y)) .* (x == 0) ...
        + (16*pi^3*exp(-t)*sin(2*pi*y) - 2*pi*exp(-t)*sin(2*pi*y)) .* (x == 1) ...
        + (2*pi*exp(-t)*sin(2*pi*x) - 16*pi^3*exp(-t)*sin(2*pi*x)) .* (y == 0) ...
        + (16*pi^3*exp(-t)*sin(2*pi*x) - 2*pi*exp(-t)*sin(2*pi*x)) .* (y == 1);
    case 'exp'
      %Vorteile für newton-version 3
      problemData.epsilon = 1;
      problemData.mobility = 'non-constant';
      problemData.mobilityCont = @(c) exp(c);
      problemData.cCont = @(t,x,y) exp(-t) * sin(2 * pi * x) .* sin(2 * pi * y);
      problemData.c0Cont = @(x,y) problemData.cCont(0,x,y);
      problemData.fCont = @(t,x,y)exp(exp(-t).*sin(x.*pi.*2.0).*sin(y.*pi.*2.0)).*(pi.^2.*exp(t.*-3.0).*sin(x.*pi.*2.0).^3.*sin(y.*pi.*2.0).^3.*8.0+pi.^4.*exp(-t).*sin(x.*pi.*2.0).*sin(y.*pi.*2.0).*3.2e1+pi.^2.*exp(-t).*sin(x.*pi.*2.0).*sin(y.*pi.*2.0).*(exp(t.*-2.0).*sin(x.*pi.*2.0).^2.*sin(y.*pi.*2.0).^2-1.0).*4.0-pi.^2.*exp(t.*-3.0).*cos(x.*pi.*2.0).^2.*sin(x.*pi.*2.0).*sin(y.*pi.*2.0).^3.*2.4e1)+exp(exp(-t).*sin(x.*pi.*2.0).*sin(y.*pi.*2.0)).*(pi.^2.*exp(t.*-3.0).*sin(x.*pi.*2.0).^3.*sin(y.*pi.*2.0).^3.*8.0+pi.^4.*exp(-t).*sin(x.*pi.*2.0).*sin(y.*pi.*2.0).*3.2e1+pi.^2.*exp(-t).*sin(x.*pi.*2.0).*sin(y.*pi.*2.0).*(exp(t.*-2.0).*sin(x.*pi.*2.0).^2.*sin(y.*pi.*2.0).^2-1.0).*4.0-pi.^2.*exp(t.*-3.0).*cos(y.*pi.*2.0).^2.*sin(x.*pi.*2.0).^3.*sin(y.*pi.*2.0).*2.4e1)-exp(-t).*sin(x.*pi.*2.0).*sin(y.*pi.*2.0)-pi.*exp(exp(-t).*sin(x.*pi.*2.0).*sin(y.*pi.*2.0)).*exp(-t).*cos(x.*pi.*2.0).*sin(y.*pi.*2.0).*(pi.^3.*exp(-t).*cos(x.*pi.*2.0).*sin(y.*pi.*2.0).*1.6e1+pi.*exp(t.*-3.0).*cos(x.*pi.*2.0).*sin(x.*pi.*2.0).^2.*sin(y.*pi.*2.0).^3.*4.0+pi.*exp(-t).*cos(x.*pi.*2.0).*sin(y.*pi.*2.0).*(exp(t.*-2.0).*sin(x.*pi.*2.0).^2.*sin(y.*pi.*2.0).^2-1.0).*2.0).*2.0-pi.*exp(exp(-t).*sin(x.*pi.*2.0).*sin(y.*pi.*2.0)).*exp(-t).*cos(y.*pi.*2.0).*sin(x.*pi.*2.0).*(pi.^3.*exp(-t).*cos(y.*pi.*2.0).*sin(x.*pi.*2.0).*1.6e1+pi.*exp(t.*-3.0).*cos(y.*pi.*2.0).*sin(x.*pi.*2.0).^3.*sin(y.*pi.*2.0).^2.*4.0+pi.*exp(-t).*cos(y.*pi.*2.0).*sin(x.*pi.*2.0).*(exp(t.*-2.0).*sin(x.*pi.*2.0).^2.*sin(y.*pi.*2.0).^2-1.0).*2.0).*2.0;
      problemData.gNContCon = @(t,x,y) ...
        - 2 * pi * exp(-t) * sin(2*pi*y) .* (x == 0) ...
        + 2 * pi * exp(-t) * sin(2*pi*y) .* (x == 1) ...
        - 2 * pi * exp(-t) * sin(2*pi*x) .* (y == 0) ...
        + 2 * pi * exp(-t) * sin(2*pi*x) .* (y == 1);
      problemData.gNContPot = @(t,x,y) ...
        + (2*pi*exp(-t)*sin(2*pi*y) - 16*pi^3*exp(-t)*sin(2*pi*y)) .* (x == 0) ...
        + (16*pi^3*exp(-t)*sin(2*pi*y) - 2*pi*exp(-t)*sin(2*pi*y)) .* (x == 1) ...
        + (2*pi*exp(-t)*sin(2*pi*x) - 16*pi^3*exp(-t)*sin(2*pi*x)) .* (y == 0) ...
        + (16*pi^3*exp(-t)*sin(2*pi*x) - 2*pi*exp(-t)*sin(2*pi*x)) .* (y == 1);
  end
elseif strcmp(problemData.scenario,'rand')
  % Randomly initialized initial state
  rng('default')
  problemData.epsilon = 2 * problemData.hmax; % Constant interface parameter
  problemData.mobility = 1; % Constant mobility parameter
  scaler = 0.99;
  problemData.c0Cont = @(x1, x2) scaler * (randi(3, size(x1,1),size(x1,2)) - 2);
  problemData.fCont = @(t,x1,x2) zeros(size(x1));
  problemData.gNContCon = @(t,x1,x2) zeros(size(x1));
  problemData.gNContPot = @(t,x1,x2) zeros(size(x1));
elseif strcmp(problemData.scenario,'ball')
  % Initial configuration is a square which transforms into a ball
%   problemData.epsilon = problemData.hmax; % Constant interface parameter
  problemData.epsilon = 1/16; % Constant interface parameter
  problemData.mobility = 1; % Constant mobility parameter
  scaler = 1.9;
  if problemData.Voronoi_tri
    problemData.c0Cont = @(x1, x2) scaler * (abs( (x1-0.75) ) < 0.2) .* (abs( (x2-0.5) ) < 0.2) - scaler/2;
  else
    problemData.c0Cont = @(x1, x2) scaler * (abs( (x1-0.5) ) < 0.25) .* (abs( (x2-0.5) ) < 0.25) - scaler/2 ...
      + 16 * scaler .* (x1 - 03/16) .* (x1 > 03/16) .* (x1 < 04/16) .* (x2 > 0.25) .* (x2 < 0.75) ...
      + 16 * scaler .* (13/16 - x1) .* (x1 > 12/16) .* (x1 < 13/16) .* (x2 > 0.25) .* (x2 < 0.75) ...
      + 16 * scaler .* (x2 - 03/16) .* (x2 > 03/16) .* (x2 < 04/16) .* (x1 > 0.25) .* (x1 < 0.75) ...
      + 16 * scaler .* (13/16 - x2) .* (x2 > 12/16) .* (x2 < 13/16) .* (x1 > 0.25) .* (x1 < 0.75);
  end
  problemData.fCont = @(t,x1,x2) zeros(size(x1));
  problemData.gNContCon = @(t,x1,x2) zeros(size(x1));
  problemData.gNContPot = @(t,x1,x2) zeros(size(x1));
elseif strcmp(problemData.scenario,'front')
  % Initial configuration is a square which transforms into a ball
  problemData.epsilon = problemData.hmax; % Constant interface parameter
  problemData.mobility = 1; % Constant mobility parameter
  mover = 2/16;
  problemData.c0Cont = @(x1, x2) -(x1 < 0.25 + mover) + (4 * (x1 - mover) - 2) .* (x1 > 0.25 + mover) .* (x1 < 0.75 + mover) + (x1 > 0.75 + mover);
  problemData.fCont = @(t,x1,x2) zeros(size(x1));
  problemData.gNContCon = @(t,x1,x2) zeros(size(x1));
  problemData.gNContPot = @(t,x1,x2) zeros(size(x1));
elseif strcmp(problemData.scenario,'merge-balls')
  % Initial configuration is two balls which transforms into a ball
  problemData.epsilon = problemData.hmax; % Constant interface parameter
%   problemData.epsilon = 1/16; % Constant interface parameter
  problemData.mobility = 1; % Constant mobility parameter
  scaler = 1.8;
  if problemData.Voronoi_tri
    problemData.c0Cont = @(x1, x2) min(0.95, scaler * ( ( (x1-7/12).^2 + (x2-sqrt(3)/6).^2 ) < 0.06 ) ...
      + scaler * ( ( (x1-11/12).^2 + (x2-sqrt(3)/3).^2 ) < 0.06 ) - scaler/2);
  else
    problemData.c0Cont = @(x1, x2) min(0.95, scaler * ( ( (x1-1/3).^2 + (x2-1/3).^2 ) < 0.0625 ) ...
      + scaler * ( ( (x1-2/3).^2 + (x2-2/3).^2 ) < 0.0625 ) - scaler/2);
  end
  problemData.fCont = @(t,x1,x2) zeros(size(x1));
  problemData.gNContCon = @(t,x1,x2) zeros(size(x1));
  problemData.gNContPot = @(t,x1,x2) zeros(size(x1));
end % problemData.scenario


% Convergence test for diffusion
% problemData.c0Cont = @(x1,x2) sin(2*pi*(x1));
% problemData.gNContCon = @(t,x1,x2) 2 * pi * cos(2*pi*(x1+t)) .* ((x1==0) - (x1==1));
% problemData.fCont = @(t,x1,x2) 4*pi*pi*sin(2*pi*(x1 + t))+2*pi*cos(2*pi*(x1+t));

% Linear test case for diffusion
% problemData.c0Cont = @(x1,x2) x1;
% problemData.gNCont = @(t,x1,x2) (x1 == 0) - (x1 == 1);
% problemData.fCont = @(t,x1,x2) 0*x1;

%% Parameter check.
assert(problemData.p >= 0 && problemData.p <= 4, 'Polynomial order must be zero to four.')
% assert(problemData.ordRK >= 1 && problemData.ordRK <= 3, 'Order of Runge Kutta must be zero to three.')
assert(problemData.hmax > 0, 'Maximum edge length must be positive.')
% assert(problemData.numSteps > 0, 'Number of time steps must be positive.')
assert(~problemData.initialLimiting || problemData.p > 0, 'Slope limiting only available for p > 0.')
assert(problemData.hmax <= problemData.epsilon, 'Grid width must be smaller than or equal to interface parameter')
assert(~(problemData.isFluxLim && ~isnumeric(problemData.mobility)), ...
  'FluxLimiting in combination with non-constant mobility is currently not supported')
%% Domain and triangulation configuration.
% Triangulate unit square using pdetool (if available or Friedrichs-Keller otherwise).
if license('checkout','PDE_Toolbox') && ~problemData.Voronoi_tri && ~problemData.Friedrichs_Keller
  problemData.generateGridData = @(hmax) domainPolygon([0 1 1 0], [0 0 1 1], hmax);
elseif problemData.Voronoi_tri
  problemData.generateGridData = @domainSquare;
else
  fprintf('PDE_Toolbox not available. Using Friedrichs-Keller triangulation.\n');
  problemData.generateGridData = @domainSquare;
end % if
problemData.generateMarkE0Tint = @(g) (g.idE0T == 0);
problemData.generateMarkE0TbdrD = @(g) false(g.numT,3);
problemData.generateMarkE0TbdrN = @(g) ~(g.markE0Tint | g.markE0TbdrD);
end % function

function psi = obstaclePotentialPsi(C)

psi = inf(size(C)); % Initialize by infinity.
tol = 1E-10;
markC = (C >= -1 - tol) & (C <= 1+tol); % Mark concentrations that are mapped to values away from infinity.
psi(markC) = max(0.5*(1-C(markC).^2), 0); 

end

function psi = obstaclePotentialPsiPrime(C)

psi = inf(size(C)); % Initialize by infinity.
psi(C<0) = -inf;
tol = 1E-10;
markC = (C >= -1 - tol) & (C <= 1+tol); % Mark concentrations that are mapped to values away from infinity.
psi(markC) = -C(markC); 

end
