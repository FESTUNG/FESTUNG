function problemData = configureProblem(problemData)
% Polynomial approximation order
p = 2;

% Configuration for shallow water solver
problemData.sweData = struct;
problemData.sweData.isCoupling = true;
problemData.sweData.p = p;

problemData.sweData = execin('swe/configureProblem', problemData.sweData);
problemData.sweData.numSteps = 3142; % Number of time steps

% Configuration for transport solver
problemData.transportData = struct;
problemData.transportData.p = p;
problemData.transportData.ordRK = problemData.sweData.schemeOrder; % as of now both models have to use the same RK method
problemData.transportData.tEnd = problemData.sweData.tEnd;
problemData.transportData.numSteps = problemData.sweData.numSteps; % as of now both models have to use the same time step
problemData.transportData.isVisGrid = false; % visualization of grid

problemData.transportData = execin('transport/configureProblem', problemData.transportData);
end % function