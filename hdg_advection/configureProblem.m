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
problemData.hmax        = 2^-4; % maximum edge length of triangle
problemData.p           =2; % local polynomial degree
problemData.ordRK       = 3; % order of Runge Kutta time stepper.
% problemData.ordRK       = min(problemData.p+1,4); % order of Runge Kutta time stepper.
problemData.numSteps    = 80; % number of time steps
problemData.tEnd        = 1; % end time

problemData.isVisGrid   = false; % visualization of grid
problemData.isVisSol    = true; % visualization of solution

% problemData.outputFrequency = max(problemData.numSteps/16,1); % no visualization of every timestep
problemData.outputFrequency = 640; % no visualization of every timestep
problemData.outputBasename  = ['output' filesep 'solution_hdg_advection']; % Basename of output files
problemData.outputTypes     = {'vtk'}; % solution output file types

%% HDG specific parameters
problemData.stab = 1.0; %stabilization parameter in mod. LF/Rusanov flux
problemData.isTrueLocalSolve = true;

%% HDG related configuration

%% Testing?
problemData.isInTesting = false;
problemData.tabRK = getDIRKtableau( problemData.ordRK );
problemData.showWaitBar = false;
problemData.showFprintfProgress = true;

%% Parameter check.
assert(problemData.p >= 0 && problemData.p <= 4, 'Polynomial order must be zero to four.')
assert(problemData.ordRK >= 1 && problemData.ordRK <= 4, 'Order of Runge Kutta must be zero to three.')
assert(problemData.hmax > 0, 'Maximum edge length must be positive.')
assert(problemData.numSteps > 0, 'Number of time steps must be positive.')
%% Coefficients and boundary data (rotating Gaussian).
problemData.getLinearAdvectionSol = @(t, X1, X2) sin( 2. * pi .* (X1 - 1. * t) ) .* sin( 2. * pi .* (X2 - 1. * t)  );


problemData.c0Cont = @(x1, x2) problemData.getLinearAdvectionSol(0, x1, x2);
problemData.fCont = @(t,x1,x2) zeros(size(x1));
problemData.cDCont = @(t,x1,x2) problemData.getLinearAdvectionSol(t, x1, x2);
problemData.gNCont = @(t,x1,x2) zeros(size(x1));

problemData.fluxCont = @( t, x1, x2, c ) evalLinearAdvectionFlux(t, x1, x2, c);

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
% problemData.generateMarkE0TbdrN = @(g) false(g.numT,3);
problemData.generateMarkE0TbdrN = @(g) generateLinearAdvectionBoundary(g);
problemData.generateMarkE0TbdrD = @(g) ~(g.markE0Tint | g.markE0TbdrN);

end % function
