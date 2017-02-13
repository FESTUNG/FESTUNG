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
%problemData.hmax        = 2^-3; % maximum edge length of triangle
problemData.hmax        = 2^-0; % maximum edge length of triangle
problemData.p           = 1; % local polynomial degree
problemData.ordRK       = min(problemData.p+1,3); % order of Runge Kutta time stepper.
problemData.numSteps    = 1; % number of time steps
problemData.tEnd        = pi/(4*16); % end time

problemData.isVisGrid   = false; % visualization of grid
problemData.isVisSol    = true; % visualization of solution

problemData.outputFrequency = 1; % no visualization of every timestep
problemData.outputBasename  = ['output' filesep 'solution_hdg_advection']; % Basename of output files
problemData.outputTypes     = {'vtk'}; % solution output file types

%% Parameter check.
assert(problemData.p >= 0 && problemData.p <= 4, 'Polynomial order must be zero to four.')
assert(problemData.ordRK >= 1 && problemData.ordRK <= 3, 'Order of Runge Kutta must be zero to three.')
assert(problemData.hmax > 0, 'Maximum edge length must be positive.')
assert(problemData.numSteps > 0, 'Number of time steps must be positive.')
%% Coefficients and boundary data (rotating Gaussian).
problemData.rgX1c = -0.2;
% problemData.rgX1c = -0.0;
problemData.rgX2c =  0.0;
problemData.rgEps =  0.0;
problemData.rgS =  0.1;
problemData.rgS2 = problemData.rgS^2;

problemData.getRGX1 = @(t, X1, X2)  X1 .* cos(4 * t) + X2 .* sin(4*t) - problemData.rgX1c;
problemData.getRGX2 = @(t, X1, X2) -X1 .* sin(4 * t) + X2 .* cos(4*t) - problemData.rgX2c;
problemData.getRGRadSq = @(t, X1, X2) problemData.getRGX1(t, X1, X2).^2 + problemData.getRGX2(t, X1, X2).^2;
problemData.getRGSol = @(t, X1, X2) (2*problemData.rgS2) / (2 * problemData.rgS2 + 4 * problemData.rgEps * t) .* exp( - ( problemData.getRGRadSq(t, X1, X2) ) ./ (2. * problemData.rgS2 + 4 * problemData.rgEps * t ) );

% problemData.getRGSol = @(t, X1, X2) zeros(size(X1));
% problemData.getRGSol = @(t, X1, X2) ones(size(X1));
problemData.getRGSol = @(t, X1, X2) X1;
problemData.u1Cont = @(t,x1,x2) zeros(size(x1));
problemData.u2Cont = @(t,x1,x2) zeros(size(x1));

problemData.c0Cont = @(x1, x2) problemData.getRGSol(0, x1, x2);
problemData.fCont = @(t,x1,x2) zeros(size(x1));
% problemData.u1Cont = @(t,x1,x2) -4.*x2;
% problemData.u2Cont = @(t,x1,x2)  4.*x1;
%problemData.cDCont = @(t,x1,x2) zeros(size(x1));
problemData.cDCont = @(t,x1,x2) problemData.getRGSol(t, x1, x2);
% problemData.cDCont = @(t,x1,x2) zeros(size(x1));
problemData.gNCont = @(t,x1,x2) zeros(size(x1));

problemData.fluxCont = @( t, x1, x2, c ) [  problemData.u1Cont(t, x1, x2) .* c; problemData.u2Cont(t, x1, x2) .* c ];


%% HDG specific parameters
problemData.stab = 1.0; %stabilization parameter in mod. LF/Rusanov flux

%% Domain and triangulation configuration.
% Triangulate unit square using pdetool (if available or Friedrichs-Keller otherwise).

problemData.generateGridData = @(hmax) domainArbitrarySquare( -0.5, 0.5, hmax );
% problemData.generateGridData = @(hmax) domainPolygon([-0.5 0.5 0.5 -0.5], [-0.5 -0.5 0.5 0.5], hmax);
% if license('checkout','PDE_Toolbox')
%   problemData.generateGridData = @(hmax) domainPolygon([-0.5 0.5 0.5 -0.5], [-0.5 -0.5 0.5 0.5], hmax);
% else
%   fprintf('PDE_Toolbox not available. Using Friedrichs-Keller triangulation.\n');
%   problemData.generateGridData = @domainSquare;
% end % if
% Specify edge ids of boundary conditions
problemData.generateMarkE0Tint = @(g) g.idE0T == 0;
problemData.generateMarkE0TbdrN = @(g) false(g.numT,3);
problemData.generateMarkE0TbdrD = @(g) ~(g.markE0Tint | g.markE0TbdrN);
end % function