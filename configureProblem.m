function problemData = configureProblem(problemData)

%% Parameters.
% Name of testcase
problemData = setdefault(problemData, 'testcase', 'convergence');

% Number of elements in x- and y-direction
problemData = setdefault(problemData, 'numElem', [32, 16]);

% Local polynomial approximation order (0 to 5)
problemData = setdefault(problemData, 'p', 2);

% Order of quadrature rule
problemData = setdefault(problemData, 'qOrd', 2*problemData.p + 1);

% Time stepping parameters
problemData = setdefault(problemData, 't0', 0);  % start time
problemData = setdefault(problemData, 'tEnd', 86.4);  % end time
problemData = setdefault(problemData, 'numSteps', 700);  % number of time steps

% Visualization settings
problemData = setdefault(problemData, 'isVisGrid', false);  % visualization of grid
problemData = setdefault(problemData, 'isVisSol', false);  % visualization of solution
problemData = setdefault(problemData, 'outputFrequency', 10); % no visualization of every timestep
problemData = setdefault(problemData, 'outputBasename', ...  % Basename of output files
                         ['output' filesep 'solution_sweVert_' problemData.testcase ]); 
problemData = setdefault(problemData, 'outputTypes', { 'vtk' });  % Type of visualization files ('vtk, 'tec')

%% Parameter check.
assert(problemData.p >= 0 && problemData.p <= 5, 'Polynomial order must be zero to five.')
assert(problemData.numSteps > 0, 'Number of time steps must be positive.')

%% Coefficients and boundary data.
[problemData, domainWidth, xi0Cont, zBotCont, idLand, idOS, idRiv, idRad] = getTestcase(problemData, problemData.testcase);
generateX = @(numElem) (0:numElem(1)) * domainWidth / numElem(1);
generateZbot = @(numElem) zBotCont(generateX(numElem));
generateXi0 = @(numElem) xi0Cont(generateX(numElem));

%% Domain and triangulation.
fn_domainRectTrap = getFunctionHandle('darcyVert/domainRectTrap');
problemData.generateGrid = @(numElem) fn_domainRectTrap(generateX(numElem), [generateZbot(numElem); generateXi0(numElem)], numElem);
problemData.generateGrid1D = @(numElem, g2D) generateGridData1D(generateX(numElem), generateXi0(numElem), numElem, g2D);

% Boundary parts (0 = int, 1 = bot, 2 = right, 3 = top, 4 = left)
checkMultipleIds = @(idE0T, ids) logical(sum(bsxfun(@eq, idE0T, reshape(ids, 1, 1, length(ids))), 3));

problemData.generateMarkE0Tint = @(g) g.idE0T == 0;
problemData.generateMarkE0TbdrBot = @(g) g.idE0T == 1;
problemData.generateMarkE0TbdrTop = @(g) g.idE0T == 3;
problemData.generateMarkE0TbdrLand = @(g) checkMultipleIds(g.idE0T, idLand);
problemData.generateMarkE0TbdrOS = @(g) checkMultipleIds(g.idE0T, idOS);
problemData.generateMarkE0TbdrRiv = @(g) checkMultipleIds(g.idE0T, idRiv);
problemData.generateMarkE0TbdrRad = @(g) checkMultipleIds(g.idE0T, idRad);

problemData.generateMarkV0T1Dint = @(g) g.idV0T == 0;
problemData.generateMarkV0T1DbdrLand = @(g) checkMultipleIds(g.idV0T, idLand);
problemData.generateMarkV0T1DbdrOS = @(g) checkMultipleIds(g.idV0T, idOS);
problemData.generateMarkV0T1DbdrRiv = @(g) checkMultipleIds(g.idV0T, idRiv);
problemData.generateMarkV0T1DbdrRad = @(g) checkMultipleIds(g.idV0T, idRad);
end % function
