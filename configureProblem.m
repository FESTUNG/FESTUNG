function problemData = configureProblem(problemData)

%% Parameters.
domainWidth = 100;  % width of computational domain

% Number of elements in x- and y-direction
problemData = setdefault(problemData, 'numElem', [16, 16]);

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
problemData = setdefault(problemData, 'isVisSol', false);  % visualization of solution
problemData = setdefault(problemData, 'outputFrequency', 10); % no visualization of every timestep
problemData = setdefault(problemData, 'outputBasename', ['output' filesep 'solution_sweVert' ]); % Basename of output files
problemData = setdefault(problemData, 'outputTypes', { 'vtk' });

%% Parameter check.
assert(problemData.p >= 0 && problemData.p <= 4, 'Polynomial order must be zero to four.')
assert(problemData.numSteps > 0, 'Number of time steps must be positive.')

%% Coefficients and boundary data.

if license('checkout', 'Symbolic_Toolbox')
  syms x z t
  
  gSym = sym('10');
  zBotSym = sym('0');
  
  DSym = { symfun(0.001, [t x z]),     symfun(0, [t x z]) ; ...
               symfun(0, [t x z]), symfun(0.001, [t x z]) };

  deltaSym = sym('0.01');
  rhoSym = sym('0.1');
         
  hSym(t,x) = deltaSym * sin(rhoSym * (t + x)) + 2;
  u1Sym(t,x,z) = sqrt(deltaSym) * z * sin(deltaSym * (t + x));
  u2Sym(t,x,z) = -0.5 * deltaSym^1.5 * z^2 * cos(deltaSym * (t + x));
  
  [problemData, h0Const, zBotConst] = analyticalData(problemData, hSym, u1Sym, u2Sym, gSym, zBotSym, DSym, domainWidth);
end % if

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
% z_bot = -0.5;
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
                     
% 
% quadratic velocity and linear height and linear diffusion
%
% z_bot = 0;
% h_0 = 1;
% h_var = 0.05;
% problemData.gConst = 10;
% % Diffusion matrix
% problemData.DCont = cell(2,2);
% problemData.DCont{1,1} = @(t,x1,x2) x1 + 1;
% problemData.DCont{1,2} = @(t,x1,x2) ones(size(x1));
% problemData.DCont{2,1} = @(t,x1,x2) ones(size(x1));
% problemData.DCont{2,2} = @(t,x1,x2) x2 + 1;
% % Analytical solution
% problemData.hCont = @(t,x1) (x1-0.5)*h_var + 1;
% problemData.u1Cont = @(t,x1,x2) x1 - x2.^2;
% problemData.u2Cont = @(t,x1,x2) -x2;
% problemData.q1Cont = @(t,x1,x2) ones(size(x1));
% problemData.q2Cont = @(t,x1,x2) -2*x2;
% % Analytical right hand side
% problemData.fhCont = @(t,x1) problemData.hCont(t,x1) + h_var * x1 - h_var * (problemData.hCont(t,x1) + z_bot).^2;
% problemData.fuCont = @(t,x1,x2) x1 + x2.^2 + problemData.gConst * h_var + 1 + 4 * x2;
                              
% 
% piecewise linear velocity and height
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
% problemData.hCont = @(t,x1) (x1-0.5)*h_var + 1 + (x1>0.5)*0.05;
% problemData.u1Cont = @(t,x1,x2) (1+x1>0.5).*x1;
% problemData.u2Cont = @(t,x1,x2) -(1+x1>0.5).*x2;
% problemData.q1Cont = @(t,x1,x2) -(1+x1>0.5).*ones(size(x1));
% problemData.q2Cont = @(t,x1,x2) zeros(size(x1));
% % Analytical right hand side
% problemData.fhCont = @(t,x1) (1+x1>0.5).*(problemData.hCont(t,x1) + h_var * x1);
% problemData.fuCont = @(t,x1,x2) (1+x1>0.5).*x1 + problemData.gConst * h_var;
             
% 
% quadratic velocity and quadratic height
%
% z_bot = 0;
% h_0 = 1;
% h_var = 0.05;
% problemData.gConst = 10;
% % Diffusion matrix
% problemData.DCont = cell(2,2);
% problemData.DCont{1,1} = @(t,x1,x2) x1 + 1;
% problemData.DCont{1,2} = @(t,x1,x2) ones(size(x1));
% problemData.DCont{2,1} = @(t,x1,x2) ones(size(x1));
% problemData.DCont{2,2} = @(t,x1,x2) x2 + 1;
% % Analytical solution
% problemData.hCont = @(t,x1) (x1.^2-0.5)*h_var + 1;
% problemData.u1Cont = @(t,x1,x2) x1 - x2.^2;
% problemData.u2Cont = @(t,x1,x2) -x2;
% problemData.q1Cont = @(t,x1,x2) ones(size(x1));
% problemData.q2Cont = @(t,x1,x2) -2*x2;
% % Analytical right hand side
% problemData.fhCont = @(t,x1) problemData.hCont(t,x1) + 2 * h_var * x1 - 2 * h_var * x1 .* (problemData.hCont(t,x1) + z_bot).^2;
% problemData.fuCont = @(t,x1,x2) x1 + x2.^2 + 2 * problemData.gConst * h_var * x1 + 1 + 4 * x2;
                       
% z_bot = 2;
% h_0 = 2;
% paramE = 0.01;
% problemData.gConst = 10;
% % Diffusion matrix
% problemData.DCont = cell(2,2);
% problemData.DCont{1,1} = @(t,x1,x2) zeros(size(x1)) + 1.e-3;
% problemData.DCont{1,2} = @(t,x1,x2) zeros(size(x1));
% problemData.DCont{2,1} = @(t,x1,x2) zeros(size(x1));
% problemData.DCont{2,2} = @(t,x1,x2) zeros(size(x1)) + 1.e-3;
% % Analytical solution
% problemData.hCont = @(t,x1) paramE * sin(paramE * (t + x1)) + h_0;
% problemData.u1Cont = @(t,x1,x2) sqrt(paramE) * (x2 - z_bot) .* sin(paramE * (t + x1));
% problemData.u2Cont = @(t,x1,x2) -0.5 * sqrt(paramE) * paramE * (x2 - z_bot).^2 .* cos(paramE * (t + x1));
% problemData.q1Cont = @(t,x1,x2) sqrt(paramE) * paramE * (x2 - z_bot) .* cos(paramE * (t + x1));
% problemData.q2Cont = @(t,x1,x2) sqrt(paramE) * sin(paramE * (t + x1));
% % Analytical right hand side
% problemData.fhCont = @(t,x1) paramE * paramE * cos(paramE * (t + x1)) - ...
%                       0.5 * sqrt(paramE) * paramE * cos(paramE * (t + x1)) .* problemData.hCont(t,x1) .* ...
%                         ( problemData.hCont(t,x1) + 2 * paramE * sin(paramE * (t + x1)) );
% problemData.fuCont = @(t,x1,x2) sqrt(paramE) * paramE * (x2 - z_bot) .* cos(paramE * (t + x1)) .* ...
%                         ( 1 + 0.5 * sqrt(paramE) * (x2 - z_bot) .* sin(paramE * (t + x1)) ) + ...
%                        problemData.gConst * paramE * paramE * cos(paramE * (t + x1)) + ...
%                        1.e-3 .* paramE * paramE * sqrt(paramE) * (x2 - z_bot) .* sin(paramE * (t + x1));
%                      
%% Domain and triangulation.
fn_domainRectTrap = getFunctionHandle('darcyVert/domainRectTrap');
problemData.generateGrid = @(numElem) fn_domainRectTrap([0, domainWidth], [zBotConst, zBotConst + h0Const], numElem);
problemData.generateGrid1D = @(numElem, g2D) generateGridData1D([0, domainWidth], zBotConst + h0Const, numElem, g2D);
end % function

function [problemData, h0Const, zBotConst] = analyticalData(problemData, hSym, u1Sym, u2Sym, gConst, zBotSym, DSym, domainWidth)
syms x z t

%% Partial derivatives of solution
dxU1Sym = diff(u1Sym, x);
dzU1Sym = diff(u1Sym, z);
dzU2Sym = diff(u2Sym, z);

%% Check continuity
assert(isequal(dxU1Sym + dzU2Sym, symfun(0, [t x z])), 'u1 and u2 do not fulfill continuity equation')

%% Compute right hand sides
fhSym = diff(hSym, t) + diff(int(u1Sym, z, zBotSym, zBotSym + hSym), x);

fuSym = diff(u1Sym,t) + u1Sym * (2 * dxU1Sym + dzU2Sym) + u2Sym * dzU1Sym + symfun(diff(gConst * hSym, x), [t x z]) - ...
        dxU1Sym * (diff(DSym{1,1}, x) + diff(DSym{2,1}, z)) - dzU1Sym * (diff(DSym{1,2}, x) + diff(DSym{2,2}, z)) - ...
        DSym{1,1} * diff(u1Sym, x, 2) - DSym{2,2} * diff(u1Sym, z, 2) - (DSym{1,2} + DSym{2,1}) * diff(dxU1Sym, z);
      
%% Create function handles
problemData.hCont = matlabFunction(hSym, 'Vars', [t x]);
problemData.u1Cont = matlabFunction(u1Sym, 'Vars', [t x z]);
problemData.u2Cont = matlabFunction(u2Sym, 'Vars', [t x z]);

problemData.fhCont = matlabFunction(fhSym, 'Vars', [t x]);
problemData.fuCont = matlabFunction(fuSym, 'Vars', [t x z]);

problemData.DCont = cellfun(@(c) matlabFunction(c, 'Vars', [t x z]), DSym, 'UniformOutput', false);
% for i = 1 : 2
%   for j = 1 : 2
%     problemData.DCont{i, j} = matlabFunction(DSym(i, j), 'Vars', [t x z]);
%   end % for j
% end % for i

%% Determine constants
problemData.gConst = double(gConst);
zBotConst = double(zBotSym);
h0Const = double(int(hSym(0, x), x, 0, domainWidth) / domainWidth);
end % function
