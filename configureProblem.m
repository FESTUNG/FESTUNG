function problemData = configureProblem(problemData)

problemData.sweData = struct;
addpath('swe');
problemData.sweData = configureProblem(problemData.sweData);
rmpath('swe');

%% Coefficients and boundary data (LeVeque's solid body rotation). % TODO
% problemData.u1Cont = @(t,x1,x2) 0.5 - x2;
% problemData.u2Cont = @(t,x1,x2) x1 - 0.5;

problemData.transportData = struct;

problemData.transportData.ordRK = problemData.sweData.schemeOrder; % as of now both models have to use the same RK method
problemData.transportData.tEnd = problemData.sweData.tEnd;
problemData.transportData.numSteps = problemData.sweData.numSteps; % as of now both models have to use the same time step

problemData.transportData.isVisGrid = false; % visualization of grid

% problemData.transportData.u1Cont = problemData.u1Cont;
% problemData.transportData.u2Cont = problemData.u2Cont;

addpath('transport');
problemData.transportData = configureProblem(problemData.transportData);
rmpath('transport');
end % function