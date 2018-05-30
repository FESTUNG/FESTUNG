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
%> @f$c:\overline{J}\times\overline{\Omega}\rightarrow\mathbb{R}@f$ 
%> of the diffusion equation
%> @f{align*}{
%> \mathbf{z}                                  &\;=\; - \nabla c   &&\text{in}~J\times\Omega\,,\\
%> \partial_t c  + \nabla\cdot (d\,\mathbf{z}) &\;=\; f            &&\text{in}~J\times\Omega\,,\\
%> c                                           &\;=\; c_\mathrm{D} &&\text{on}~J\times{\partial\Omega}_{\mathrm{D}}\,,\\
%> \vec{z}\cdot\vec{\nu}                       &\;=\; g_\mathrm{N} &&\text{on}~J\times{\partial\Omega}_\mathrm{N}\,,\\
%> c                                           &\;=\; c^0          &&\text{on}~\{0\}\times\Omega\,.
%> @f}
%> The vector-valued quantity @f$\mathbf{z}@f$ was introduced as auxiliary unknown.  
%> The coefficients @f$d:J\times\Omega\rightarrow\mathbb{R}^+@f$ and 
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
problemData = setdefault(problemData, 'testcase', 'showcase');  % select testcase
problemData = setdefault(problemData, 'hmax', 2^-3);            % maximum edge length of triangle
problemData = setdefault(problemData, 'p', 2);                  % local polynomial degree
problemData = setdefault(problemData, 'tEnd', pi);              % end time
problemData = setdefault(problemData, 'numSteps', 20);          % number of time steps
problemData = setdefault(problemData, 'isVisGrid', false);      % visualization of grid
problemData = setdefault(problemData, 'isVisSol', true);        % visualization of solution
problemData = setdefault(problemData, 'eta', 1);                % penalty parameter (eta>0)
problemData = setdefault(problemData, 'outputBasename', ['output' filesep 'solution_diffusion']); % basename of output files
problemData = setdefault(problemData, 'outputTypes', {'vtk', 'tec'}); % type of output files
%% Parameter check.
assert(problemData.p >= 0 && problemData.p <= 4, 'Polynomial order must be zero to four.')
assert(problemData.hmax > 0, 'Maximum edge length must be positive.' )
assert(problemData.numSteps > 0, 'Number of time steps must be positive.')
%% Coefficients and boundary data.
switch problemData.testcase
  case 'showcase'
    problemData.isStationary = false;
    
    problemData.c0Cont = @(x1,x2) sin(x1) .* cos(x2);
    problemData.dCont  = @(t,x1,x2) (x1 < 3/4 & x1 > 1/4 & x2 < 3/4 & x2 > 1/4) + 0.01;
    problemData.fCont  = @(t,x1,x2) 0.1 * t * (x1 == x1);
    problemData.cDCont = @(t,x1,x2) sin(2 * pi * x2 + t);
    problemData.gNCont = @(t,x1,x2) x2;
    
  case 'convergence'
    problemData.isStationary = true;
    
    problemData.cCont = @(t,x1,x2) cos(7 * x1) .* cos(7 * x2);
    problemData.dCont = @(t,x1,x2) exp(x1 + x2);
    problemData.fCont = @(t,x1,x2) 98 * exp(x1 + x2) .* cos(7 * x1) .* cos(7 * x2) ...
      + 7 * exp(x1 + x2) .* (sin(7 * x1) .* cos(7 * x2) + cos(7 * x1) .* sin(7 * x2));
    problemData.fCont = @(t,x1,x2) 98 * cos(7 * x1) .* cos(7 * x2);
    problemData.c0Cont = @(x1,x2) problemData.cCont(0,x1,x2);
    problemData.cDCont = problemData.cCont;
    problemData.gNCont = @(t,x1,x2) 7 * cos(7 * x1) .* sin(7 * x2) .* ((x2 > 0.5) - (x2 < 0.5));
    
  otherwise
    error('Invalid testcase "%s".', problemData.testcase);
end % siwtch
%% Domain and triangulation configuration.
% Select domain and triangulation
problemData.generateGridData = @domainSquare;
% Specify edge ids of boundary conditions
problemData.generateMarkE0Tint = @(g) g.idE0T == 0;
problemData.generateMarkE0TbdrN = @(g) g.idE0T == 1 | g.idE0T == 3;
problemData.generateMarkE0TbdrD = @(g) ~(g.markE0Tint | g.markE0TbdrN);
end % function