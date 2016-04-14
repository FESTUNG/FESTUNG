function problemData = configureProblem(problemData)
%% Parameters.
problemData.hmax            = 2^-3;       % maximum edge length of triangle
problemData.p               = 2;          % local polynomial degree
problemData.tEnd            = pi;         % end time
problemData.numSteps        = 20;         % number of time steps
problemData.isVisGrid       = false;      % visualization of grid
problemData.isVisSol        = false;      % visualization of solution
problemData.eta             = 1;          % penalty parameter (eta>0)
problemData.outputBasename  = 'solution'; % Basename of output files
problemData.outputTypes     = cellstr(['vtk';'tec']);
%% Parameter check.
assert(problemData.p >= 0 && problemData.p <= 4, 'Polynomial order must be zero to four.')
assert(problemData.hmax > 0, 'Maximum edge length must be positive.' )
assert(problemData.numSteps > 0, 'Number of time steps must be positive.')
%% Coefficients and boundary data.
problemData.c0Cont = @(x1,x2) sin(x1).*cos(x2);
problemData.dCont  = @(t,x1,x2) (x1<3/4&x1>1/4&x2<3/4&x2>1/4) + 0.01;
problemData.fCont  = @(t,x1,x2) 0.1*t*(x1==x1);
problemData.cDCont = @(t,x1,x2) sin(2*pi*x2 + t);
problemData.gNCont = @(t,x1,x2) x2;
end % function