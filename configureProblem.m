%===============================================================================
%> @file
%>
%> @brief Sets all problem parameters and fills the problemData-struct with
%>        all basic configuration options.
%===============================================================================
%>
%> @brief Sets all problem parameters and fills the problemData-struct with
%>        all basic configuration options.
%>
%> This routine is called before any other function for the problem.
%> It defines all problem parameters and should be the only file users have
%> to adapt to their needs.
%> 
%> For an outline of the model, see @link swe_2dv @endlink.
%>
%> The actual testcases are defined in swe_2dv/getTestcase.m
%>
%> A detailed description can be found in @ref RRAFK2018.
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
%> @copyright 2014-2018 Balthasar Reuter, Florian Frank, Vadym Aizinger
%>
%> @author Balthasar Reuter, 2018
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
% Name of testcase
problemData = setdefault(problemData, 'testcase', 'showcase');

% Number of elements in x- and y-direction
problemData = setdefault(problemData, 'numElem', [32, 8]);

% Local polynomial approximation order (0 to 5)
problemData = setdefault(problemData, 'p', 1);

% Order of quadrature rule
problemData = setdefault(problemData, 'qOrd', 2*problemData.p + 1);
problemData = setdefault(problemData, 'qOrdMax', problemData.qOrd);

% Time stepping parameters
problemData = setdefault(problemData, 't0', 0);  % start time
problemData = setdefault(problemData, 'tEnd', 600);  % end time
problemData = setdefault(problemData, 'numSteps', ceil(problemData.tEnd/0.01));  % number of time steps

% Order of Runge-Kutta method (1 - explicit Euler, 2/3 - multi-stage RK)
problemData = setdefault(problemData, 'ordRK', 1);
% problemData = setdefault(problemData, 'ordRK', min(problemData.p+1, 3));

% Visualization settings
problemData = setdefault(problemData, 'isVisGrid', false);  % visualization of grid
problemData = setdefault(problemData, 'isVisSol', true);  % visualization of solution
problemData = setdefault(problemData, 'outputFrequency', 200); % no visualization of every timestep
problemData = setdefault(problemData, 'outputBasename', ...  % Basename of output files
                         ['output' filesep problemData.problemName '_' problemData.testcase ]); 
problemData = setdefault(problemData, 'outputTypes', { 'vtk' });  % Type of visualization files ('vtk, 'tec')

% Coupling settings
problemData = setdefault(problemData, 'isCoupling', false);  % Enable coupling
problemData = setdefault(problemData, 'isJumpCoupling', true);  % Include jump penalty term on coupling interface

% Reduce order of q, u1 to be p/2
problemData = setdefault(problemData, 'isReducedOrder', false);

% Bottom friction parameterization
problemData = setdefault(problemData, 'bottomFriction', 'none'); % 'none' (use uDCont), 'linear', 'quadratic'

%% Parameter check.
assert(problemData.p >= 0 && problemData.p <= 5, 'Polynomial order must be zero to five.')
assert(problemData.ordRK >= 1 && problemData.ordRK <= 3, 'Order of Runge-Kutta method must be one to three.')
assert(problemData.numSteps > 0, 'Number of time steps must be positive.')

%% Coefficients and boundary data.
problemData = execin([problemData.problemName filesep 'getTestcase'], problemData, problemData.testcase);

%% Domain and triangulation.
domainWidth = problemData.domainWidth;
generateX = @(numElem) (0:numElem(1)) * domainWidth / numElem(1);

zBotCont = problemData.zBotCont;
xi0Cont = problemData.xi0Cont;

problemData.generateGrid = @(numElem) domainRectTrap(generateX(numElem), [zBotCont(generateX(numElem)); xi0Cont(generateX(numElem))], numElem);
problemData.generateGrid1D = @(numElem, g2D) generateGridData1D(generateX(numElem), xi0Cont(generateX(numElem)), numElem, g2D);

% Boundary parts (0 = int, 1 = bot, 2 = right, 3 = top, 4 = left)
checkMultipleIds = @(idE0T, ids) logical(sum(bsxfun(@eq, idE0T, reshape(ids, 1, 1, length(ids))), 3));

idBdrRiem = problemData.idBdrRiem;
idBdrH = problemData.idBdrH;
idBdrU = problemData.idBdrU;
idBdrQ = problemData.idBdrQ;

problemData.generateMarkE0Tint = @(g) g.idE0T == 0;
problemData.generateMarkE0TbdrRiem = @(g) checkMultipleIds(g.idE0T, idBdrRiem);
problemData.generateMarkE0TbdrCoupling = @(g) g.idE0T == 1;
problemData.generateMarkE0TbdrBot = @(g) g.idE0T == 1;
problemData.generateMarkE0TbdrTop = @(g) g.idE0T == 3;
problemData.generateMarkE0TbdrH = @(g) checkMultipleIds(g.idE0T, idBdrH);
problemData.generateMarkE0TbdrU = @(g) checkMultipleIds(g.idE0T, idBdrU);
problemData.generateMarkE0TbdrQ = @(g) checkMultipleIds(g.idE0T, idBdrQ);
end % function
