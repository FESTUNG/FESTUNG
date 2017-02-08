function problemData = configureProblem(problemData)

%% Parameters.
domainWidth = 100;  % width of computational domain
problemData.numElem = 32*[1,1];%[24, 24];  % number of elements per direction
problemData.p = 1; % local polynomial degree
problemData.qOrd = 2*problemData.p + 1; % order of quadrature formula
problemData.t0 = 0; % start time
problemData.tEnd = 0.1; % end time
problemData.numSteps = 200; % number of time steps
problemData.isVisGrid = false; % visualization of grid
problemData.isVisSol = false; % visualization of solution
problemData.eta = 1; % penalty parameter (eta>0)
problemData.outputFrequency = 10; % no visualization of every timestep
problemData.outputBasename = ['output' filesep 'solution_sweVert' ]; % Basename of output files
problemData.outputTypes = { 'vtk' };

%% Parameter check.
assert(problemData.p >= 0 && problemData.p <= 4, 'Polynomial order must be zero to four.')
assert(problemData.numSteps > 0, 'Number of time steps must be positive.')

%% Coefficients and boundary data.
%
% Constant solution
%
% z_bot = 0;
% h_0 = 1;
% problemData.gConst = 10;
% % Diffusion matrix
% problemData.DCont = cell(2,2);
% problemData.DCont{1,1} = @(t,x1,x2) zeros(size(x1)) + 1.e-3;
% problemData.DCont{1,2} = @(t,x1,x2) zeros(size(x1));
% problemData.DCont{2,1} = @(t,x1,x2) zeros(size(x1));
% problemData.DCont{2,2} = @(t,x1,x2) zeros(size(x1)) + 1.e-3;
% % Analytical solution
% problemData.hCont = @(t,x1) ones(size(x1));
% problemData.u1Cont = @(t,x1,x2) ones(size(x1));
% problemData.u2Cont = @(t,x1,x2) zeros(size(x1));
% problemData.q1Cont = @(t,x1,x2) zeros(size(x1));
% problemData.q2Cont = @(t,x1,x2) zeros(size(x1));
% % Analytical right hand side
% problemData.fhCont = @(t,x1) zeros(size(x1));
% problemData.fuCont = @(t,x1,x2) zeros(size(x1));
           
% 
% Linear height
%
% z_bot = 0;
% h_0 = 1;
% h_var = 0.05;
% problemData.gConst = 10;
% % Diffusion matrix
% problemData.DCont = cell(2,2);
% problemData.DCont{1,1} = @(t,x1,x2) zeros(size(x1)) + 1.e-3;
% problemData.DCont{1,2} = @(t,x1,x2) zeros(size(x1));
% problemData.DCont{2,1} = @(t,x1,x2) zeros(size(x1));
% problemData.DCont{2,2} = @(t,x1,x2) zeros(size(x1)) + 1.e-3;
% % Analytical solution
% problemData.hCont = @(t,x1) (x1-0.5)*h_var + 1;
% problemData.u1Cont = @(t,x1,x2) ones(size(x1));
% problemData.u2Cont = @(t,x1,x2) zeros(size(x1));
% problemData.q1Cont = @(t,x1,x2) zeros(size(x1));
% problemData.q2Cont = @(t,x1,x2) zeros(size(x1));
% % Analytical right hand side
% problemData.fhCont = @(t,x1) ones(size(x1))*h_var;
% problemData.fuCont = @(t,x1,x2)  problemData.gConst*h_var * ones(size(x1));
%    
% 
% z-Linear velocity
%
% z_bot = 0;
% h_0 = 1;
% problemData.gConst = 10;
% % Diffusion matrix
% problemData.DCont = cell(2,2);
% problemData.DCont{1,1} = @(t,x1,x2) zeros(size(x1)) + 1.e-3;
% problemData.DCont{1,2} = @(t,x1,x2) zeros(size(x1));
% problemData.DCont{2,1} = @(t,x1,x2) zeros(size(x1));
% problemData.DCont{2,2} = @(t,x1,x2) zeros(size(x1)) + 1.e-3;
% % Analytical solution
% problemData.hCont = @(t,x1) ones(size(x1));
% problemData.u1Cont = @(t,x1,x2) x2;
% problemData.u2Cont = @(t,x1,x2) zeros(size(x1));
% problemData.q1Cont = @(t,x1,x2) zeros(size(x1));
% problemData.q2Cont = @(t,x1,x2) -ones(size(x1));
% % Analytical right hand side
% problemData.fhCont = @(t,x1) zeros(size(x1));
% problemData.fuCont = @(t,x1,x2) zeros(size(x1));
                 
% 
% x-Linear velocity
%
% z_bot = 0;
% h_0 = 1;
% problemData.gConst = 10;
% % Diffusion matrix
% problemData.DCont = cell(2,2);
% problemData.DCont{1,1} = @(t,x1,x2) zeros(size(x1)) + 1.e-3;
% problemData.DCont{1,2} = @(t,x1,x2) zeros(size(x1));
% problemData.DCont{2,1} = @(t,x1,x2) zeros(size(x1));
% problemData.DCont{2,2} = @(t,x1,x2) zeros(size(x1)) + 1.e-3;
% % Analytical solution
% problemData.hCont = @(t,x1) ones(size(x1));
% problemData.u1Cont = @(t,x1,x2) x1;
% problemData.u2Cont = @(t,x1,x2) -x2;
% problemData.q1Cont = @(t,x1,x2) -ones(size(x1));
% problemData.q2Cont = @(t,x1,x2) zeros(size(x1));
% % Analytical right hand side
% problemData.fhCont = @(t,x1) ones(size(x1));
% problemData.fuCont = @(t,x1,x2) x1;
                       
% 
% linear velocity and height
%
% z_bot = 0;
% h_0 = 1;
% h_var = 0.05;
% problemData.gConst = 10;
% % Diffusion matrix
% problemData.DCont = cell(2,2);
% problemData.DCont{1,1} = @(t,x1,x2) zeros(size(x1)) + 1.e-3;
% problemData.DCont{1,2} = @(t,x1,x2) zeros(size(x1));
% problemData.DCont{2,1} = @(t,x1,x2) zeros(size(x1));
% problemData.DCont{2,2} = @(t,x1,x2) zeros(size(x1)) + 1.e-3;
% % Analytical solution
% problemData.hCont = @(t,x1) (x1-0.5)*h_var + 1;
% problemData.u1Cont = @(t,x1,x2) x1;
% problemData.u2Cont = @(t,x1,x2) -x2;
% problemData.q1Cont = @(t,x1,x2) -ones(size(x1));
% problemData.q2Cont = @(t,x1,x2) zeros(size(x1));
% % Analytical right hand side
% problemData.fhCont = @(t,x1) problemData.hCont(t,x1) + h_var * x1;
% problemData.fuCont = @(t,x1,x2) x1 + problemData.gConst * h_var;
   
% 
% Quadratic height
%
% z_bot = 0;
% h_0 = 1;
% h_var = 0.1;
% problemData.gConst = 10;
% % Diffusion matrix
% problemData.DCont = cell(2,2);
% problemData.DCont{1,1} = @(t,x1,x2) zeros(size(x1)) + 1.e-3;
% problemData.DCont{1,2} = @(t,x1,x2) zeros(size(x1));
% problemData.DCont{2,1} = @(t,x1,x2) zeros(size(x1));
% problemData.DCont{2,2} = @(t,x1,x2) zeros(size(x1)) + 1.e-3;
% % Analytical solution
% problemData.hCont = @(t,x1) (x1.^2-0.5)*h_var + 1;
% problemData.u1Cont = @(t,x1,x2) ones(size(x1));
% problemData.u2Cont = @(t,x1,x2) zeros(size(x1));
% problemData.q1Cont = @(t,x1,x2) zeros(size(x1));
% problemData.q2Cont = @(t,x1,x2) zeros(size(x1));
% % Analytical right hand side
% problemData.fhCont = @(t,x1) 2*h_var*x1;
% problemData.fuCont = @(t,x1,x2)  2*problemData.gConst*h_var * x1;
   

% 
% quadratic velocity
%
% z_bot = 0;
% h_0 = 1;
% h_var = 0.05;
% problemData.gConst = 10;
% % Diffusion matrix
% problemData.DCont = cell(2,2);
% problemData.DCont{1,1} = @(t,x1,x2) ones(size(x1));
% problemData.DCont{1,2} = @(t,x1,x2) zeros(size(x1));
% problemData.DCont{2,1} = @(t,x1,x2) zeros(size(x1));
% problemData.DCont{2,2} = @(t,x1,x2) ones(size(x1));
% % Analytical solution
% problemData.hCont = @(t,x1) ones(size(x1));
% problemData.u1Cont = @(t,x1,x2) x1 - x2.^2;
% problemData.u2Cont = @(t,x1,x2) -x2;
% problemData.q1Cont = @(t,x1,x2) ones(size(x1));
% problemData.q2Cont = @(t,x1,x2) -2*x2;
% % Analytical right hand side
% problemData.fhCont = @(t,x1) problemData.hCont(t,x1);
% problemData.fuCont = @(t,x1,x2) x1 + x2.^2 + 2;
                  
% 
% quadratic velocity and linear height
%
% z_bot = 0;
% h_0 = 1;
% h_var = 0.05;
% problemData.gConst = 10;
% % Diffusion matrix
% problemData.DCont = cell(2,2);
% problemData.DCont{1,1} = @(t,x1,x2) ones(size(x1));
% problemData.DCont{1,2} = @(t,x1,x2) zeros(size(x1));
% problemData.DCont{2,1} = @(t,x1,x2) zeros(size(x1));
% problemData.DCont{2,2} = @(t,x1,x2) ones(size(x1));
% % Analytical solution
% problemData.hCont = @(t,x1) (x1-0.5)*h_var + 1;
% problemData.u1Cont = @(t,x1,x2) x1 - x2.^2;
% problemData.u2Cont = @(t,x1,x2) -x2;
% problemData.q1Cont = @(t,x1,x2) ones(size(x1));
% problemData.q2Cont = @(t,x1,x2) -2*x2;
% % Analytical right hand side
% problemData.fhCont = @(t,x1) problemData.hCont(t,x1) + h_var * x1 - h_var * (problemData.hCont(t,x1) + z_bot).^2;
% problemData.fuCont = @(t,x1,x2) x1 + x2.^2 + 2 + problemData.gConst * h_var;
                       
z_bot = 2;
h_0 = 2;
paramE = 0.01;
problemData.gConst = 10;
% Diffusion matrix
problemData.DCont = cell(2,2);
problemData.DCont{1,1} = @(t,x1,x2) zeros(size(x1)) + 1.e-3;
problemData.DCont{1,2} = @(t,x1,x2) zeros(size(x1));
problemData.DCont{2,1} = @(t,x1,x2) zeros(size(x1));
problemData.DCont{2,2} = @(t,x1,x2) zeros(size(x1)) + 1.e-3;
% Analytical solution
problemData.hCont = @(t,x1) paramE * sin(paramE * (t + x1)) + h_0;
problemData.u1Cont = @(t,x1,x2) sqrt(paramE) * (x2 - z_bot) .* sin(paramE * (t + x1));
problemData.u2Cont = @(t,x1,x2) -0.5 * sqrt(paramE) * paramE * (x2 - z_bot).^2 .* cos(paramE * (t + x1));
problemData.q1Cont = @(t,x1,x2) sqrt(paramE) * paramE * (x2 - z_bot) .* cos(paramE * (t + x1));
problemData.q2Cont = @(t,x1,x2) sqrt(paramE) * sin(paramE * (t + x1));
% Analytical right hand side
problemData.fhCont = @(t,x1) paramE * paramE * cos(paramE * (t + x1)) - ...
                      0.5 * sqrt(paramE) * paramE * cos(paramE * (t + x1)) .* problemData.hCont(t,x1) .* ...
                        ( problemData.hCont(t,x1) + 2 * paramE * sin(paramE * (t + x1)) );
problemData.fuCont = @(t,x1,x2) sqrt(paramE) * paramE * (x2 - z_bot) .* cos(paramE * (t + x1)) .* ...
                        ( 1 + 0.5 * sqrt(paramE) * (x2 - z_bot) .* sin(paramE * (t + x1)) ) + ...
                       problemData.gConst * paramE * paramE * cos(paramE * (t + x1)) + ...
                       1.e-3 .* paramE * paramE * sqrt(paramE) * (x2 - z_bot) .* sin(paramE * (t + x1));
                     
%% Domain and triangulation.
fn_domainRectTrap = getFunctionHandle('darcyVert/domainRectTrap');
problemData.generateGrid = @(numElem) fn_domainRectTrap([0, domainWidth], [z_bot, z_bot + h_0], numElem);
problemData.generateGrid1D = @(numElem, g2D) generateGridData1D([0, domainWidth], z_bot + h_0, numElem, g2D);
end % function
