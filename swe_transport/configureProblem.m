function problemData = configureProblem(problemData)

problemData.sweData = struct; % TODO
problemData.transportData = struct;

%% Parameters.
problemData.hmax        = 2^-6; % maximum edge length of triangle
problemData.p           = 1; % local polynomial degree
problemData.ordRK       = min(problemData.p+1,3); % order of Runge Kutta time stepper.
problemData.tEnd        = 100/3142*2*pi; % end time 
problemData.numSteps    = 100/3142*3142; % number of time steps

%% Coefficients and boundary data (LeVeque's solid body rotation).
problemData.u1Cont = @(t,x1,x2) 0.5 - x2;
problemData.u2Cont = @(t,x1,x2) x1 - 0.5;

problemData.transportData.ordRK = problemData.ordRK; % as of now both models have to use the same RK method
problemData.transportData.tEnd = problemData.tEnd;
problemData.transportData.numSteps = problemData.numSteps; % as of now both models have to use the same time step

problemData.transportData.isVisGrid = false; % visualization of grid

problemData.transportData.u1Cont = problemData.u1Cont;
problemData.transportData.u2Cont = problemData.u2Cont;

% problemData.sweData = ...

%% Add problem to search path
addpath('transport');
problemData.transportData = configureProblem(problemData.transportData);
rmpath('transport');
end % function