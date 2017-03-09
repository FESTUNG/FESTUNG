function problemData = configureProblem(problemData)

%% Parameters.
domainWidth = 100;  % width of computational domain

% Number of elements in x- and y-direction
problemData = setdefault(problemData, 'numElem', [16, 4]);

% Local polynomial approximation order (0 to 4)
problemData = setdefault(problemData, 'p', 1);

% Order of quadrature rule
problemData = setdefault(problemData, 'qOrd', 2*problemData.p + 1);

% Time stepping parameters
problemData = setdefault(problemData, 't0', 0);  % start time
problemData = setdefault(problemData, 'tEnd', 0.1);  % end time
problemData = setdefault(problemData, 'numSteps', 10);  % number of time steps

% Visualization settings
problemData = setdefault(problemData, 'isVisGrid', false);  % visualization of grid
problemData = setdefault(problemData, 'isVisSol', true);  % visualization of solution
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
else
  error('Symbolic Toolbox required to derive problem formulation!')
end % if

%% Domain and triangulation.
fn_domainRectTrap = getFunctionHandle('darcyVert/domainRectTrap');
problemData.generateGrid = @(numElem) fn_domainRectTrap([0, domainWidth], [zBotConst, zBotConst + h0Const], numElem);
problemData.generateGrid1D = @(numElem, g2D) generateGridData1D([0, domainWidth], zBotConst + h0Const, numElem, g2D);

% Boundary parts (0 = int, 1 = bot, 2 = right, 3 = top, 4 = left)
idLand = -1; idOS = -1; idRiv = -1; idRad = -1;

problemData.generateMarkE0Tint = @(g) g.idE0T == 0;
problemData.generateMarkE0TbdrBot = @(g) g.idE0T == 1;
problemData.generateMarkE0TbdrTop = @(g) g.idE0T == 3;
problemData.generateMarkE0TbdrLand = @(g) g.idE0T == idLand;
problemData.generateMarkE0TbdrOS = @(g) g.idE0T == idOS;
problemData.generateMarkE0TbdrRiv = @(g) g.idE0T == idRiv;
problemData.generateMarkE0TbdrRad = @(g) g.idE0T == idRad;

problemData.generateMarkV0T1Dint = @(g) g.idV0T == 0;
problemData.generateMarkV0T1DbdrLand = @(g) g.idV0T == idLand;
problemData.generateMarkV0T1DbdrOS = @(g) g.idV0T == idOS;
problemData.generateMarkV0T1DbdrRiv = @(g) g.idV0T == idRiv;
problemData.generateMarkV0T1DbdrRad = @(g) g.idV0T == idRad;
end % function

function [problemData, h0Const, zBotConst] = analyticalData(problemData, hSym, u1Sym, u2Sym, gConst, zBotSym, DSym, domainWidth)
syms x z t

%% Partial derivatives of solution
dxU1Sym = diff(u1Sym, x);
dzU1Sym = diff(u1Sym, z);
dzU2Sym = diff(u2Sym, z);

%% Check continuity
assert(isequal(dxU1Sym + dzU2Sym, symfun(0, [t x z])), 'u1 and u2 do not fulfill continuity equation')

%% Compute boundary conditions
qDSym = -sign(x - 0.5 * domainWidth) * (DSym{1,1} * dxU1Sym + DSym{1,2} * dzU1Sym);
depthIntU1Sym = int(u1Sym, z);
depthIntU1Cont = matlabFunction(depthIntU1Sym, 'Vars', [t x z]);

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

problemData.hDCont = problemData.hCont;
problemData.u1DCont = problemData.u1Cont;
problemData.u2DCont = problemData.u2Cont;
problemData.qDCont = matlabFunction(qDSym, 'Vars', [t x z]);
problemData.uhDCont = @(t, x, xi, zb) problemData.hDCont(t, x) .* (depthIntU1Cont(t, x, xi) - depthIntU1Cont(t, x, zb)) ./ (xi - zb);

%% Determine constants
problemData.gConst = double(gConst);
zBotConst = double(zBotSym);
h0Const = double(int(hSym(0, x), x, 0, domainWidth) / domainWidth);
end % function
