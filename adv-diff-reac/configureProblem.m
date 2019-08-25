% Fills the problemData-struct with all basic configuration options.
% Problem parameters are to be modified inside this routine.

%===============================================================================
%> @file
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
%> @f$u:\overline{J}\times\overline{\Omega}\rightarrow\mathbb{R}@f$ 
%> of the advection-diffusion-reaction equation
%> @f{align*}{
%> \partial_t u + \nabla \cdot (\vec{v}\,u - d\,\nabla u) + r\,u  &\;=\; 0            && \text{in}~J\times\Omega\,, \\
%>                                                              u &\;=\; u_\mathrm{D} && \text{on}~J\times\partial\Omega_\mathrm{D}\,, \\
%>                                    -d \,\nabla u \cdot \vec{n} &\;=\; g_\mathrm{N} && \text{on}~J\times\partial\Omega_\mathrm{N}\,, \\
%>                                                              u &\;=\; u^0          && \text{on}~\{t_0\}\times\Omega\
%> @f}
%> with given initial $u^0:\Omega\to\mathbb{R}^+_0@f$ and boundary data 
%> @f$u_\mathrm{D}: J\times\partial\Omega_\mathrm{D} \to \mathbb{R}^+_0@f$,
%> @f$g_\mathrm{N}: J\times\partial\Omega_\mathrm{N} \to \mathbb{R}@f$.
%> 
%> A detailed description can be found in @ref RHRAFK2019.
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
%> @copyright 2014-2019 Balthasar Reuter, Florian Frank, Vadym Aizinger
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
problemData = setdefault(problemData, 'hmax', 2^-3);            % maximum edge length of triangle
problemData = setdefault(problemData, 'p', 1);                  % local polynomial degree
problemData = setdefault(problemData, 't0', 0);                 % start time
problemData = setdefault(problemData, 'tEnd', .05);               % end time
problemData = setdefault(problemData, 'numSteps', 5e3);         % number of time steps
problemData = setdefault(problemData, 'isVisGrid', false);      % visualization of grid
problemData = setdefault(problemData, 'isVisSol', true);        % visualization of solution
problemData = setdefault(problemData, 'symparam', 1);           % symmetrization parameter (theta)
problemData = setdefault(problemData, 'penparam', 10);           % penalty parameter (eta>0)
problemData = setdefault(problemData, 'isIP', true);           % use interior penalty dG instead of LDG
problemData = setdefault(problemData, 'deltaAdvReac', 0);       % time stepping parameter for advection (0, 1)
problemData = setdefault(problemData, 'deltaDiff', 0);          % time stepping parameter for diffusion (0, 1)
problemData = setdefault(problemData, 'outputBasename', ['output' filesep 'adv-diff-reac']); % basename of output files
problemData = setdefault(problemData, 'outputTypes', {'vtk'});  % type of output files
problemData = setdefault(problemData, 'outputFrequency', 1e3);  % output frequency
%% Coefficients and boundary data.
uCont = @(t,x1,x2) cos(7*x1) .* cos(7*x2) + 3*t;
problemData.uCont = uCont;
problemData.u0Cont = @(x1,x2) uCont(0,x1,x2);
problemData.uDCont = uCont;
problemData.gNCont = @(t,x1,x2) zeros(size(x1));
problemData.rCont = @(t,x1,x2) sin(5*x1) .* sin(5*x2);
problemData.dCont = @(t,x1,x2) exp(0.5 * (x1 + x2));
problemData.v1Cont = @(t,x1,x2) -sin(3*x1) .* sin(3*x2) + 2;
problemData.v2Cont = @(t,x1,x2) -cos(3*x1) .* cos(3*x2) + 2;
problemData.fCont = @(t,x1,x2) 7 * ( (sin(3*x1) .* sin(3*x2) + 2 - 0.5 * exp(0.5 * (x1+x2))) .* sin(7*x1) .* cos(7*x2) + ...
                                     (cos(3*x1) .* cos(3*x2) + 2 - 0.5 * exp(0.5 * (x1+x2))) .* cos(7*x1) .* cos(7*x2) ) - ...
                               98 * exp(0.5 * (x1+x2)) .* cos(7*x1) .* cos(7*x2) + ...
                               sin(5*x1) .* sin(5*x2) .* (cos(7*x1) .* cos(7*x2) + 3*t);
%% Domain and triangulation configuration.
% Select domain and triangulation
problemData.generateGridData = @domainSquare;
% Specify edge ids of boundary conditions
problemData.generateMarkE0Tint = @(g) g.idE0T == 0;
problemData.generateMarkE0TbdrN = @(g) g.idE0T == 4 | g.idE0T == 3; %false(g.numT,3);
problemData.generateMarkE0TbdrD = @(g) ~(g.markE0Tint | g.markE0TbdrN);
end % function