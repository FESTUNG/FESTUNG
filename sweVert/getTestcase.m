function [problemData, domainWidth, h0Const, zBotConst, idLand, idOS, idRiv, idRad] = getTestcase(problemData, name)
switch name
  case 'constant'
    domainWidth = 1;
    idLand = -1; idOS = -1; idRiv = -1; idRad = -1;
    
    problemData.gConst = 10;
    h0Const = 1;
    zBotConst = 0;
    
    problemData.hCont = @(t,x) ones(size(x));
    problemData.u1Cont = @(t,x,z) ones(size(x));
    problemData.u2Cont = @(t,x,z) zeros(size(x));
    
    problemData.fhCont = @(t,x) zeros(size(x));
    problemData.fuCont = @(t,x,z) zeros(size(x));
    
    problemData.DCont = { @(t,x,z) ones(size(x)), @(t,x,z) zeros(size(x)); 
                          @(t,x,z) zeros(size(x)), @(t,x,z) ones(size(x)) };
    
    problemData.hDCont = problemData.hCont;
    problemData.u1DCont = problemData.u1Cont;
    problemData.u2DCont = problemData.u2Cont;
    problemData.qDCont = @(t,x,z) zeros(size(x));
    problemData.uhDCont = @(t,x) ones(size(x));

  case 'linear h'
    domainWidth = 1;
    idLand = -1; idOS = -1; idRiv = -1; idRad = -1;
    
    problemData.gConst = 10;
    h0Const = 1;
    zBotConst = 0;
    hVar = 0.05;
    
    problemData.hCont = @(t,x) 1 + (x - domainWidth/2) * hVar;
    problemData.u1Cont = @(t,x,z) ones(size(x));
    problemData.u2Cont = @(t,x,z) zeros(size(x));
    
    problemData.fhCont = @(t,x) hVar * ones(size(x));
    problemData.fuCont = @(t,x,z) problemData.gConst * hVar * ones(size(x));
    
    problemData.DCont = { @(t,x,z) ones(size(x)), @(t,x,z) zeros(size(x)); 
                          @(t,x,z) zeros(size(x)), @(t,x,z) ones(size(x)) };
    
    problemData.hDCont = problemData.hCont;
    problemData.u1DCont = problemData.u1Cont;
    problemData.u2DCont = problemData.u2Cont;
    problemData.qDCont = @(t,x,z) zeros(size(x));
    problemData.uhDCont = @(t,x) problemData.hCont(t,x);

  case 'z-linear u'
    domainWidth = 1;
    idLand = -1; idOS = -1; idRiv = -1; idRad = -1;
    
    problemData.gConst = 10;
    h0Const = 1;
    zBotConst = 0;
    
    problemData.hCont = @(t,x) ones(size(x));
    problemData.u1Cont = @(t,x,z) z;
    problemData.u2Cont = @(t,x,z) zeros(size(x));
    
    problemData.fhCont = @(t,x) zeros(size(x));
    problemData.fuCont = @(t,x,z) zeros(size(x));
    
    problemData.DCont = { @(t,x,z) ones(size(x)), @(t,x,z) zeros(size(x)); 
                          @(t,x,z) zeros(size(x)), @(t,x,z) ones(size(x)) };
    
    problemData.hDCont = problemData.hCont;
    problemData.u1DCont = problemData.u1Cont;
    problemData.u2DCont = problemData.u2Cont;
    problemData.qDCont = @(t,x,z) zeros(size(x));
    problemData.uhDCont = @(t,x) 0.5 * (problemData.hCont(t,x).^2 - zBotConst.^2);
  
  case 'x-linear u'
    domainWidth = 1;
    idLand = -1; idOS = -1; idRiv = -1; idRad = -1;
    
    problemData.gConst = 10;
    h0Const = 1;
    zBotConst = 0;
    
    problemData.hCont = @(t,x) ones(size(x));
    problemData.u1Cont = @(t,x,z) x;
    problemData.u2Cont = @(t,x,z) -z;
    
    problemData.fhCont = @(t,x) ones(size(x));
    problemData.fuCont = @(t,x,z) x;
    
    problemData.DCont = { @(t,x,z) ones(size(x)), @(t,x,z) zeros(size(x)); 
                          @(t,x,z) zeros(size(x)), @(t,x,z) ones(size(x)) };
    
    problemData.hDCont = problemData.hCont;
    problemData.u1DCont = problemData.u1Cont;
    problemData.u2DCont = problemData.u2Cont;
    problemData.qDCont = @(t,x,z) zeros(size(x));
    problemData.uhDCont = @(t,x) problemData.hCont(t,x) .* problemData.u1Cont(t,x,0);

  case 'linear'
    domainWidth = 1;
    idLand = -1; idOS = -1; idRiv = -1; idRad = -1;
    
    problemData.gConst = 10;
    h0Const = 1;
    zBotConst = 0;
    hVar = 0.05;
    
    problemData.hCont = @(t,x) 1 + (x - domainWidth/2) * hVar;
    problemData.u1Cont = @(t,x,z) x;
    problemData.u2Cont = @(t,x,z) -z;
    
    problemData.fhCont = @(t,x) problemData.hCont(t,x) + hVar * x;
    problemData.fuCont = @(t,x,z) x + z.^2 + problemData.gConst * hVar + 2;
    
    problemData.DCont = { @(t,x,z) ones(size(x)), @(t,x,z) zeros(size(x)); 
                          @(t,x,z) zeros(size(x)), @(t,x,z) ones(size(x)) };
    
    problemData.hDCont = problemData.hCont;
    problemData.u1DCont = problemData.u1Cont;
    problemData.u2DCont = problemData.u2Cont;
    problemData.qDCont = @(t,x,z) zeros(size(x));
    problemData.uhDCont = @(t,x) problemData.hCont(t,x) .* problemData.u1Cont(t,x,0);

  case 'convergence'
    if license('checkout', 'Symbolic_Toolbox')
      syms x z t

      domainWidth = 100;
      idLand = -1; idOS = -1; idRiv = 2; idRad = 4;

      gSym = sym('10');
      zBotSym = sym('0');

      deltaSym = sym('0.01');
      rhoSym = sym('0.1');

      hSym(t,x) = deltaSym * sin(rhoSym * (t + x)) + 2;
      u1Sym(t,x,z) = sqrt(deltaSym) * z * sin(deltaSym * (t + x));
      u2Sym(t,x,z) = -0.5 * deltaSym^1.5 * z^2 * cos(deltaSym * (t + x));
      DSym = cellfun(@(c) symfun(c, [t x z]), {0.001, 0; 0, 0.001}, 'UniformOutput', false);

      [problemData, h0Const, zBotConst] = analyticalData(problemData, hSym, u1Sym, u2Sym, gSym, zBotSym, DSym, domainWidth);
    else
      error('Symbolic Toolbox required to derive problem formulation!')
    end % if
end % switch
end % function


function [problemData, h0Const, zBotConst] = analyticalData(problemData, hSym, u1Sym, u2Sym, gConst, zBotSym, DSym, domainWidth)
syms x z t

%% Partial derivatives of solution
dxU1Sym = diff(u1Sym, x);
dzU1Sym = diff(u1Sym, z);
dzU2Sym = diff(u2Sym, z);

%% Check continuity
assert(isequal(dxU1Sym + dzU2Sym, symfun(0, [t x z])), 'u1 and u2 do not fulfill continuity equation')

%% Depth integrated velocity
depthIntU1Sym = int(u1Sym, z, zBotSym, zBotSym + hSym);

%% Compute boundary conditions
qDSym = -sign(x - 0.5 * domainWidth) * (DSym{1,1} * dxU1Sym + DSym{1,2} * dzU1Sym);

%% Compute right hand sides
fhSym = diff(hSym, t) + diff(depthIntU1Sym, x);

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
problemData.uhDCont = matlabFunction(depthIntU1Sym, 'Vars', [t x]);

%% Determine constants
problemData.gConst = double(gConst);
zBotConst = double(zBotSym);
h0Const = double(int(hSym(0, x), x, 0, domainWidth) / domainWidth);
end % function
