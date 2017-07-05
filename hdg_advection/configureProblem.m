% Fills the problemData-struct with all basic configuration options.
% Problem parameters are to be modified inside this routine.

%===============================================================================
%> @file advection/configureProblem.m
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
%> @f$c:\overline{J}\times\overline{\Omega}\rightarrow\mathbb{R}@f$ 
%> of the advection equation
%> @f{align*}{
%> \partial_t c  + \nabla\cdot (\mathbf{u}\,c) &\;=\; f            &&\text{in}~J\times\Omega\,,\\
%> c                                           &\;=\; c_\mathrm{D} &&\text{on}~J\times{\partial\Omega}_{\mathrm{in}}\,,\\
%> c                                           &\;=\; c^0          &&\text{on}~\{0\}\times\Omega\,.
%> @f}
%> The velocity @f$\mathbf{u}:J\times\Omega\rightarrow\mathbb{R}^2@f$ and 
%> right hand side@f$f:J\times\Omega\rightarrow \mathbb{R}@f$
%> may vary in time and space. 
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
problemData = setdefault(problemData, 'hmax', 2^-4); % maximum edge length of triangle
problemData = setdefault(problemData, 'p', 3); % local polynomial degree
% problemData = setdefault(problemData, 'ordRK', 2); % order of Runge-Kutta time stepper.
problemData.ordRK       = min(problemData.p+1,4); % order of Runge Kutta time stepper.
% problemData.ordRK       = 1; % order of Runge Kutta time stepper.
problemData = setdefault(problemData, 'numSteps', 160); % number of time steps
% problemData = setdefault(problemData, 'tEnd', 1); % end time
problemData = setdefault(problemData, 'tEnd', 2*pi); % end time

problemData = setdefault(problemData, 'isVisGrid', false);
problemData = setdefault(problemData, 'isVisSol', true);

problemData = setdefault(problemData, 'outputFrequency', 10);
problemData = setdefault(problemData, 'outputBasename', ['output' filesep 'solution_hdg_advection']);
problemData = setdefault(problemData, 'outputTypes', {'vtk'});

%% HDG specific parameters
problemData = setdefault(problemData, 'stab', 1.0);
problemData = setdefault(problemData, 'isTrueLocalSolve', true);

%% HDG related configuration

%% Testing?
% problemData.isInTesting = false;
problemData.tabRK = getDIRKtableau( problemData.ordRK );

problemData = setdefault(problemData, 'isInTesting', false);
problemData = setdefault(problemData, 'showWaitBar', false);
problemData = setdefault(problemData, 'showFprintfProgress', false);

%% Parameter check.
assert(problemData.p >= 0 && problemData.p <= 4, 'Polynomial order must be zero to four.')
assert(problemData.ordRK >= 1 && problemData.ordRK <= 4, 'Order of Runge Kutta must be zero to three.')
assert(problemData.hmax > 0, 'Maximum edge length must be positive.')
assert(problemData.numSteps > 0, 'Number of time steps must be positive.')
%% Coefficients and boundary data (rotating Gaussian).
% problemData.getLinearAdvectionSol = @(t, X1, X2) sin( 2. * pi .* (X1 - 1. * t) ) .* sin( 2. * pi .* (X2 - 1. * t)  );
% problemData.c0Cont = @(x1, x2) problemData.getLinearAdvectionSol(0, x1, x2);
% problemData.fCont = @(t,x1,x2) zeros(size(x1));
% problemData.cDCont = @(t,x1,x2) problemData.getLinearAdvectionSol(t, x1, x2);
% problemData.gNCont = @(t,x1,x2) zeros(size(x1));
% problemData.fluxCont = @( t, x1, x2, c ) evalLinearAdvectionFlux(t, x1, x2, c);

%% Coefficients and boundary data (LeVeque's solid body rotation).
G = @(x1, x2, x1_0, x2_0) (1/0.15) * sqrt((x1-x1_0).^2 + (x2-x2_0).^2);
problemData.c0Cont = @(x1, x2) ((x1 - 0.5).^2 + (x2 - 0.75).^2 <= 0.0225 & (x1 <= 0.475 | x1 >= 0.525 | x2 >= 0.85)) + ...
                    (1-G(x1, x2, 0.5, 0.25)) .* ((x1 - 0.5).^2 + (x2 - 0.25).^2 <= 0.0225) + ...
                    0.25*(1+cos(pi*G(x1, x2, 0.25, 0.5))).*((x1 - 0.25).^2 + (x2 - 0.5).^2 <= 0.0225);
problemData.fCont = @(t,x1,x2) zeros(size(x1));
% problemData.u1Cont = @(t,x1,x2) 0.5 - x2;
% problemData.u2Cont = @(t,x1,x2) x1 - 0.5;
problemData.fluxCont = @( t, x1, x2, c ) evalLeVequeFlux(t, x1, x2, c);
problemData.cDCont = @(t,x1,x2) zeros(size(x1));
problemData.gNCont = @(t,x1,x2) zeros(size(x1));

problemData.getLinearAdvectionSol = @(t, X1, X2) problemData.c0Cont(X1, X2);
% 
% 
% problemData.getLinearAdvectionSol = @(t, x1, x2) cos(7*x1).*cos(7*x2);
% problemData.cCont  = @(t,x1,x2) cos(7*x1).*cos(7*x2);
% problemData.c0Cont = @(x1,x2) problemData.cCont(0,x1,x2);
% problemData.cDCont = @(t,x1,x2) problemData.cCont(0,x1,x2);
% problemData.fCont  = @(x1,x2) -7*sin(7*x1).*cos(7*x2).*exp((x1+x2)/2) ...
%                     -7*cos(7*x1).*sin(7*x2).*exp((x1-x2)/2) ...
%                     + 0.5*cos(7*x1).*cos(7*x2).*exp((x1+x2)/2) ...
%                     - 0.5*cos(7*x1).*cos(7*x2).*exp((x1-x2)/2);
% problemData.fluxCont = @(t,  x1, x2, c ) evalSteadyFlux(0, x1, x2, c);


% problemData.c0Cont = @(x1, x2) x1;
% problemData.fCont = @(x1,x2) ones(size(x1)); %source
% problemData.cDCont = @(t, x1,x2) x1;
% problemData.fluxCont = @( t, x1, x2, c ) evalLinearAdvectionFlux(t, x1, x2, c);


%% Domain and triangulation configuration.
% Triangulate unit square using pdetool (if available or Friedrichs-Keller otherwise).

% problemData.generateGridData = @(hmax) domainArbitrarySquare( -0.5, 0.5, hmax );

problemData.generateGridData = @(hmax) domainArbitrarySquare( 0.0, 1.0, hmax );

% problemData.generateGridData = @(hmax) domainPolygon([-0.5 0.5 0.5 -0.5], [-0.5 -0.5 0.5 0.5], hmax);
% if license('checkout','PDE_Toolbox')
%   problemData.generateGridData = @(hmax) domainPolygon([-0.5 0.5 0.5 -0.5], [-0.5 -0.5 0.5 0.5], hmax);
% else
%   fprintf('PDE_Toolbox not available. Using Friedrichs-Keller triangulation.\n');
%   problemData.generateGridData = @domainSquare;
% end % if
% Specify edge ids of boundary conditions
% problemData.generateMarkE0Tint = @(g) g.idE0T == 0;
% problemData.generateMarkE0TbdrN = @(g) generateLinearAdvectionBoundary(g);
% problemData.generateMarkE0TbdrD = @(g) ~(g.markE0Tint | g.markE0TbdrN);

% Specify edge ids of boundary conditions
problemData.generateMarkE0Tint = @(g) g.idE0T == 0;
% problemData.generateMarkE0TbdrN = @(g) false(g.numT,3);
problemData.generateMarkE0TbdrN = @(g) generateLeVequeBoundary( g );
problemData.generateMarkE0TbdrD = @(g) ~(g.markE0Tint | g.markE0TbdrN);

end % function
