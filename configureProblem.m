function problemData = configureProblem(problemData)

%% Parameters.
% Name of testcase
problemData = setdefault(problemData, 'testcase', 'coupling');

% Number of elements in x- and y-direction
problemData = setdefault(problemData, 'numElem', [8, 4]);

% Local polynomial approximation order (0 to 5)
problemData = setdefault(problemData, 'p', 1);

% Order of quadrature rule
problemData = setdefault(problemData, 'qOrd', 2*problemData.p + 1);

% Time stepping parameters
problemData = setdefault(problemData, 't0', 0);  % start time
problemData = setdefault(problemData, 'tEnd', 86.4);  % end time
problemData = setdefault(problemData, 'numSteps', ceil(problemData.tEnd/0.005));  % number of time steps

% Order of Runge-Kutta method
problemData = setdefault(problemData, 'ordRK', 1);%min(problemData.p+1, 3));

% Visualization settings
problemData = setdefault(problemData, 'isVisGrid', false);  % visualization of grid
problemData = setdefault(problemData, 'isVisSol', true);  % visualization of solution
problemData = setdefault(problemData, 'outputFrequency', 100); % no visualization of every timestep
problemData = setdefault(problemData, 'outputBasename', ...  % Basename of output files
                         ['output' filesep problemData.problemName '_' problemData.testcase ]); 
problemData = setdefault(problemData, 'outputTypes', { 'vtk' });  % Type of visualization files ('vtk, 'tec')

% ID of coupling boundary
problemData = setdefault(problemData, 'isCoupling', false);

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
