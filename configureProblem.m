% Fills the problemData-struct with all basic configuration options.
% Problem parameters are to be modified inside this routine.

%===============================================================================
%> @file advection-diffusion/configureProblem.m
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
%> of the diffusion equation
%> @f{align*}{
%> \mathbf{z}                                        &\;=\; - \nabla c   &&\text{in}~J\times\Omega\,,\\
%> \partial_t c  + \nabla\cdot (d\,\mathbf{z} + G c) &\;=\; f            &&\text{in}~J\times\Omega\,,\\
%> c                                                 &\;=\; c_\mathrm{D} &&\text{on}~J\times{\partial\Omega}_{\mathrm{D}}\,,\\
%> \vec{z}\cdot\vec{\nu}                             &\;=\; g_\mathrm{N} &&\text{on}~J\times{\partial\Omega}_\mathrm{N}\,,\\
%> c                                                 &\;=\; c^0          &&\text{on}~\{0\}\times\Omega\,.
%> @f}
%> The vector-valued quantity @f$\mathbf{z}@f$ was introduced as auxiliary unknown.  
%> The coefficients @f$d:J\times\Omega\rightarrow\mathbb{R}^+@f$,
%> @f$G: J\times\Omega\rightarrow\mathbb{R}^2@f$, and 
%> @f$f:J\times\Omega\rightarrow \mathbb{R}@f$ may vary in time and space. 
%> A detailed description can be found in @ref FRAK2015.
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

% Local polynomial approximation order (0 to 4)
problemData = setdefault(problemData, 'p', 1);

% Penalty parameter (eta > 0)
problemData.eta             = 1;

% Order of Runge-Kutta method
problemData = setdefault(problemData, 'ordRK', min(problemData.p+1,3));

% Time stepping parameters
problemData = setdefault(problemData, 'numSteps', 1E3);  % number of time steps
problemData = setdefault(problemData, 'tEnd', 2*pi);% (problemData.numSteps/3142)*2*pi);  % end time

% Time stepping method and reduction of scheme
problemData = setdefault(problemData, 'implicitEuler', true); % switch between implicit and explicit Euler
problemData = setdefault(problemData, 'reduceAnsatz', true); % reduce ansatz functions for flux

% Slope limiting settings (only for explicit Euler)
problemData = setdefault(problemData, 'isSlopeLim', false); % enable/disable slope limiting
problemData = setdefault(problemData, 'typeSlopeLim', 'hierarch_vert'); % Type of slope limiter (linear, hierarch_vert, strict)

% Visualization settings
problemData = setdefault(problemData, 'isVisGrid', false);  % visualization of grid
problemData = setdefault(problemData, 'isVisSol', true);  % visualization of solution
problemData = setdefault(problemData, 'outputFrequency', 100); % no visualization of every timestep
problemData = setdefault(problemData, 'outputBasename', ...
                         ['output' filesep 'advection_diffusion_' problemData.typeSlopeLim]); 
problemData = setdefault(problemData, 'outputTypes', { 'vtk' });  % Type of visualization files ('vtk, 'tec')

%% Parameter check.
assert(problemData.p >= 0 && problemData.p <= 4, 'Polynomial order must be zero to four.')
assert(problemData.ordRK >= 1 && problemData.ordRK <= 3, 'Order of Runge Kutta must be zero to three.')
assert(problemData.hmax > 0, 'Maximum edge length must be positive.')
assert(problemData.numSteps > 0, 'Number of time steps must be positive.')
assert(~problemData.isSlopeLim || problemData.p > 0, 'Slope limiting only available for p > 0.')
%% Coefficients and boundary data (LeVeque's solid body rotation).
G = @(x1, x2, x1_0, x2_0) (1/0.15) * sqrt((x1-x1_0).^2 + (x2-x2_0).^2);
problemData.dCont  = @(t,x1,x2) 0 * (x1<3/4&x1>1/4&x2<3/4&x2>1/4) + 0.0001;
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
problemData.generateMarkE0TbdrN = @(g) g.idE0T ~= 0;
problemData.generateMarkE0TbdrD = @(g) ~(g.markE0Tint | g.markE0TbdrN);
end % function