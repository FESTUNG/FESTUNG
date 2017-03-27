function problemData = configureProblem(problemData)

%% Parameters.
% Name of testcase
problemData = setdefault(problemData, 'testcase', 'convergence2');

% Number of elements in x- and y-direction
problemData = setdefault(problemData, 'numElem', [16, 8]);

% Local polynomial approximation order (0 to 5)
problemData = setdefault(problemData, 'p', 1);

% Order of quadrature rule
problemData = setdefault(problemData, 'qOrd', 2*problemData.p + 1);

% Time stepping parameters
problemData = setdefault(problemData, 't0', 0);  % start time
problemData = setdefault(problemData, 'tEnd', 0.1);  % end time
problemData = setdefault(problemData, 'numSteps', 10);  % number of time steps

% Order of Runge-Kutta methode
problemData = setdefault(problemData, 'ordRK', 1);%max(problemData.p + 1, 3));

% Visualization settings
problemData = setdefault(problemData, 'isVisGrid', false);  % visualization of grid
problemData = setdefault(problemData, 'isVisSol', false);  % visualization of solution
problemData = setdefault(problemData, 'outputFrequency', 100); % no visualization of every timestep
problemData = setdefault(problemData, 'outputBasename', ...  % Basename of output files
                         ['output' filesep 'solution_sweVert_' problemData.testcase ]); 
problemData = setdefault(problemData, 'outputTypes', { 'vtk' });  % Type of visualization files ('vtk, 'tec')

%% Parameter check.
assert(problemData.p >= 0 && problemData.p <= 5, 'Polynomial order must be zero to five.')
assert(problemData.numSteps > 0, 'Number of time steps must be positive.')

%% Coefficients and boundary data.
[problemData, domainWidth, xi0Cont, zBotCont, idLand, idOS, idRiv, idRad] = getTestcase(problemData, problemData.testcase);
problemData.h0Cont = @(x1) problemData.hCont(problemData.t0, x1);
problemData.u10Cont = @(x1,x2) problemData.u1Cont(problemData.t0, x1, x2);
generateX = @(numElem) (0:numElem(1)) * domainWidth / numElem(1);
generateZbot = @(numElem) zBotCont(generateX(numElem));
generateXi0 = @(numElem) xi0Cont(generateX(numElem));

%% Domain and triangulation.
problemData.generateGrid = @(numElem) domainRectTrap(generateX(numElem), [generateZbot(numElem); generateXi0(numElem)], numElem);
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
