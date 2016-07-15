function problemData = configureProblem(problemData)

problemData.sweData = struct;
h = getFunctionHandle('swe/configureProblem');
problemData.sweData = h(problemData.sweData);
problemData.sweData.isCoupling = true;

problemData.transportData = struct;

problemData.transportData.ordRK = problemData.sweData.schemeOrder; % as of now both models have to use the same RK method
problemData.transportData.tEnd = problemData.sweData.tEnd;
problemData.transportData.numSteps = problemData.sweData.numSteps; % as of now both models have to use the same time step

problemData.transportData.isVisGrid = false; % visualization of grid

h = getFunctionHandle('transport/configureProblem');
problemData.transportData = h(problemData.transportData);
end % function