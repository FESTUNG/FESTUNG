function problemData = configureProblem(problemData)

%% Parameters.
% Name of testcase
problemData = setdefault(problemData, 'testcase', 'linear h');

% Number of elements in x- and y-direction
problemData = setdefault(problemData, 'numElem', [8, 8]);

% Local polynomial approximation order (0 to 4)
problemData = setdefault(problemData, 'p', 1);

% Order of quadrature rule
problemData = setdefault(problemData, 'qOrd', 2*problemData.p + 1);

% Time stepping parameters
problemData = setdefault(problemData, 't0', 0);  % start time
problemData = setdefault(problemData, 'tEnd', 0.1);  % end time
problemData = setdefault(problemData, 'numSteps', 1);  % number of time steps

% Visualization settings
problemData = setdefault(problemData, 'isVisGrid', false);  % visualization of grid
problemData = setdefault(problemData, 'isVisSol', true);  % visualization of solution
problemData = setdefault(problemData, 'outputFrequency', 10); % no visualization of every timestep
problemData = setdefault(problemData, 'outputBasename', ...
                ['output' filesep 'solution_sweVert_' problemData.testcase ]); % Basename of output files
problemData = setdefault(problemData, 'outputTypes', { 'vtk' });  % Type of visualization files ('vtk, 'tec')

%% Parameter check.
assert(problemData.p >= 0 && problemData.p <= 4, 'Polynomial order must be zero to four.')
assert(problemData.numSteps > 0, 'Number of time steps must be positive.')

%% Coefficients and boundary data.

[problemData, domainWidth, h0Const, zBotConst, idLand, idOS, idRiv, idRad] = getTestcase(problemData, problemData.testcase);

%% Domain and triangulation.
fn_domainRectTrap = getFunctionHandle('darcyVert/domainRectTrap');
problemData.generateGrid = @(numElem) fn_domainRectTrap([0, domainWidth], [zBotConst, zBotConst + h0Const], numElem);
problemData.generateGrid1D = @(numElem, g2D) generateGridData1D([0, domainWidth], zBotConst + h0Const, numElem, g2D);

% Boundary parts (0 = int, 1 = bot, 2 = right, 3 = top, 4 = left)
problemData.generateMarkE0Tint = @(g) g.idE0T == 0;
problemData.generateMarkE0TbdrBot = @(g) g.idE0T == 1;
problemData.generateMarkE0TbdrTop = @(g) g.idE0T == 3;
problemData.generateMarkE0TbdrLand = @(g) g.idE0T == idLand;
problemData.generateMarkE0TbdrOS = @(g) g.idE0T == idOS;
problemData.generateMarkE0TbdrRiv = @(g) g.idE0T == idRiv;
problemData.generateMarkE0TbdrRad = @(g) g.idE0T == idRad;

problemData.generateMarkV0T1Dint = @(g) g.idV0T == 0;
problemData.generateMarkV0T1DbdrLand = @(g) g.idV0T == idLand;
problemData.generateMarkV0T1DbdrOS = @(g) g.idV0T == idOS;
problemData.generateMarkV0T1DbdrRiv = @(g) g.idV0T == idRiv;
problemData.generateMarkV0T1DbdrRad = @(g) g.idV0T == idRad;
end % function
