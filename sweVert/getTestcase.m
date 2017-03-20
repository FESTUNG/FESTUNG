function [problemData, domainWidth, xi0Cont, zBotCont, idLand, idOS, idRiv, idRad] = getTestcase(problemData, name)
switch name
  case 'constant'
    domainWidth = 1;
    idLand = -1; idOS = -1; idRiv = -1; idRad = -1;
    
    problemData.gConst = 10;
    xi0Cont = @(x) ones(size(x));
    zBotCont = @(x) zeros(size(x));
    
    problemData.hCont = @(t,x) ones(size(x));
    problemData.u1Cont = @(t,x,z) ones(size(x));
    problemData.u2Cont = @(t,x,z) zeros(size(x));
    
    problemData.fhCont = @(t,x) zeros(size(x));
    problemData.fuCont = @(t,x,z) zeros(size(x));
    
    problemData.DCont = { @(t,x,z) ones(size(x)), @(t,x,z) zeros(size(x)); ...
                          @(t,x,z) zeros(size(x)), @(t,x,z) ones(size(x)) };
    
    problemData.hDCont = problemData.hCont;
    problemData.u1DCont = problemData.u1Cont;
    problemData.u2DCont = problemData.u2Cont;
    problemData.q1DCont = @(t,x,z) zeros(size(x));
    problemData.q2DCont = @(t,x,z) zeros(size(x));
    problemData.uhDCont = @(t,x) ones(size(x));

  case 'linear h'
    domainWidth = 1;
    idLand = -1; idOS = -1; idRiv = -1; idRad = -1;
    
    problemData.gConst = 10;
    xi0Cont = @(x) ones(size(x));
    zBotCont = @(x) zeros(size(x));
    hVar = 0.05;
    
    problemData.hCont = @(t,x) 1 + (x - domainWidth/2) * hVar;
    problemData.u1Cont = @(t,x,z) ones(size(x));
    problemData.u2Cont = @(t,x,z) zeros(size(x));
    
    problemData.fhCont = @(t,x) hVar * ones(size(x));
    problemData.fuCont = @(t,x,z) problemData.gConst * hVar * ones(size(x));
    
    problemData.DCont = { @(t,x,z) ones(size(x)), @(t,x,z) zeros(size(x)); ...
                          @(t,x,z) zeros(size(x)), @(t,x,z) ones(size(x)) };
    
    problemData.hDCont = problemData.hCont;
    problemData.u1DCont = problemData.u1Cont;
    problemData.u2DCont = problemData.u2Cont;
    problemData.q1DCont = @(t,x,z) zeros(size(x));
    problemData.q2DCont = @(t,x,z) zeros(size(x));
    problemData.uhDCont = @(t,x) problemData.hCont(t,x);
    
  case 'linear zb'
    domainWidth = 1;
    idLand = -1; idOS = -1; idRiv = -1; idRad = -1;
    
    problemData.gConst = 10;
    hVar = 0.05;
    xi0Cont = @(x) zeros(size(x));
    zBotCont = @(x) -(1 + (x - domainWidth/2) * hVar);
    
    problemData.hCont = @(t,x) 1 + (x - domainWidth/2) * hVar;
    problemData.u1Cont = @(t,x,z) ones(size(x));
    problemData.u2Cont = @(t,x,z) zeros(size(x));
    
    problemData.fhCont = @(t,x) hVar * ones(size(x));
    problemData.fuCont = @(t,x,z) problemData.gConst * hVar * ones(size(x));
    
    problemData.DCont = { @(t,x,z) ones(size(x)), @(t,x,z) zeros(size(x)); ...
                          @(t,x,z) zeros(size(x)), @(t,x,z) ones(size(x)) };
    
    problemData.hDCont = problemData.hCont;
    problemData.u1DCont = problemData.u1Cont;
    problemData.u2DCont = problemData.u2Cont;
    problemData.q1DCont = @(t,x,z) zeros(size(x));
    problemData.q2DCont = @(t,x,z) zeros(size(x));
    problemData.uhDCont = @(t,x) problemData.hCont(t,x);

  case 'z-linear u'
    domainWidth = 1;
    idLand = -1; idOS = -1; idRiv = -1; idRad = -1;
    
    problemData.gConst = 10;
    xi0Cont = @(x) ones(size(x));
    zBotCont = @(x) zeros(size(x));
    
    problemData.hCont = @(t,x) ones(size(x));
    problemData.u1Cont = @(t,x,z) z;
    problemData.u2Cont = @(t,x,z) zeros(size(x));
    
    problemData.fhCont = @(t,x) zeros(size(x));
    problemData.fuCont = @(t,x,z) zeros(size(x));
    
    problemData.DCont = { @(t,x,z) ones(size(x)), @(t,x,z) zeros(size(x)); ...
                          @(t,x,z) zeros(size(x)), @(t,x,z) ones(size(x)) };
    
    problemData.hDCont = problemData.hCont;
    problemData.u1DCont = problemData.u1Cont;
    problemData.u2DCont = problemData.u2Cont;
    problemData.q1DCont = @(t,x,z) zeros(size(x));
    problemData.q2DCont = @(t,x,z) -ones(size(x));
    problemData.uhDCont = @(t,x) 0.5 * ((problemData.hCont(t,x) + zBotConst).^2 - zBotConst.^2);
  
  case 'x-linear u'
    domainWidth = 1;
    idLand = -1; idOS = -1; idRiv = -1; idRad = -1;
    
    problemData.gConst = 10;
    xi0Cont = @(x) ones(size(x));
    zBotCont = @(x) zeros(size(x));
    
    problemData.hCont = @(t,x) ones(size(x));
    problemData.u1Cont = @(t,x,z) x;
    problemData.u2Cont = @(t,x,z) -z;
    
    problemData.fhCont = @(t,x) problemData.hCont(t,x);
    problemData.fuCont = @(t,x,z) x;
    
    problemData.DCont = { @(t,x,z) ones(size(x)), @(t,x,z) zeros(size(x)); ...
                          @(t,x,z) zeros(size(x)), @(t,x,z) ones(size(x)) };
    
    problemData.hDCont = problemData.hCont;
    problemData.u1DCont = problemData.u1Cont;
    problemData.u2DCont = problemData.u2Cont;
    problemData.q1DCont = @(t,x,z) -ones(size(x));
    problemData.q2DCont = @(t,x,z) zeros(size(x));
    problemData.uhDCont = @(t,x) problemData.hCont(t,x) .* x;
    
  case 'linear h, duo-linear u'
    domainWidth = 1;
    idLand = -1; idOS = -1; idRiv = -1; idRad = -1;
    
    problemData.gConst = 10;
    xi0Cont = @(x) ones(size(x));
    zBotCont = @(x) zeros(size(x));
    hVar = 0.05;
    
    problemData.hCont = @(t,x) 1 + (x - domainWidth/2) * hVar;
    problemData.u1Cont = @(t,x,z) x + z;
    problemData.u2Cont = @(t,x,z) x - z;
    
    problemData.fhCont = @(t,x) problemData.hCont(t,x) + hVar * x + hVar * problemData.hCont(t,x);
    problemData.fuCont = @(t,x,z) 2 * x + problemData.gConst * hVar;
    
    problemData.DCont = { @(t,x,z) ones(size(x)), @(t,x,z) zeros(size(x)); ...
                          @(t,x,z) zeros(size(x)), @(t,x,z) ones(size(x)) };
    
    problemData.hDCont = problemData.hCont;
    problemData.u1DCont = problemData.u1Cont;
    problemData.u2DCont = problemData.u2Cont;
    problemData.q1DCont = @(t,x,z) -ones(size(x));
    problemData.q2DCont = @(t,x,z) -ones(size(x));
    problemData.uhDCont = @(t,x) problemData.hCont(t,x) .* x + problemData.hCont(t,x).^2 / 2;

  case 'linear'
    domainWidth = 1;
    idLand = -1; idOS = -1; idRiv = -1; idRad = -1;
    
    problemData.gConst = 10;
    xi0Cont = @(x) ones(size(x));
    zBotCont = @(x) zeros(size(x));
    hVar = 0.05;
    
    problemData.hCont = @(t,x) 1 + (x - domainWidth/2) * hVar;
    problemData.u1Cont = @(t,x,z) x;
    problemData.u2Cont = @(t,x,z) -z;
    
    problemData.fhCont = @(t,x) problemData.hCont(t,x) + hVar * x;
    problemData.fuCont = @(t,x,z) x + problemData.gConst * hVar;
    
    problemData.DCont = { @(t,x,z) ones(size(x)), @(t,x,z) zeros(size(x)); ...
                          @(t,x,z) zeros(size(x)), @(t,x,z) ones(size(x)) };
    
    problemData.hDCont = problemData.hCont;
    problemData.u1DCont = problemData.u1Cont;
    problemData.u2DCont = problemData.u2Cont;
    problemData.q1DCont = @(t,x,z) -ones(size(x));
    problemData.q2DCont = @(t,x,z) zeros(size(x));
    problemData.uhDCont = @(t,x) problemData.hCont(t,x) .* x;
    
  case 'linear h, quadratic u'
    domainWidth = 1;
    idLand = -1; idOS = -1; idRiv = -1; idRad = -1;
    
    problemData.gConst = 10;
    xi0Cont = @(x) ones(size(x));
    zBotCont = @(x) zeros(size(x));
    hVar = 0.05;
    
    problemData.hCont = @(t,x) 1 + (x - domainWidth/2) * hVar;
    problemData.u1Cont = @(t,x,z) x - z.^2;
    problemData.u2Cont = @(t,x,z) -z;
    
    problemData.fhCont = @(t,x) problemData.hCont(t,x) + hVar * x - problemData.hCont(t,x).^2 * hVar;
    problemData.fuCont = @(t,x,z) x + z.^2 + problemData.gConst * hVar + 2;
    
    problemData.DCont = { @(t,x,z) ones(size(x)), @(t,x,z) zeros(size(x)); ...
                          @(t,x,z) zeros(size(x)), @(t,x,z) ones(size(x)) };
    
    problemData.hDCont = problemData.hCont;
    problemData.u1DCont = problemData.u1Cont;
    problemData.u2DCont = problemData.u2Cont;
    problemData.q1DCont = @(t,x,z) -ones(size(x));
    problemData.q2DCont = @(t,x,z) 2 * z;
    problemData.uhDCont = @(t,x) problemData.hCont(t,x) .* x - problemData.hCont(t,x).^3 / 3;

  case 'quadratic h'
    domainWidth = 1;
    idLand = -1; idOS = -1; idRiv = -1; idRad = -1;
    
    problemData.gConst = 10;
    xi0Cont = @(x) ones(size(x));
    zBotCont = @(x) zeros(size(x));
    hVar = 0.05;
    
    problemData.hCont = @(t,x) 1 + (x.^2 - domainWidth/2) * hVar;
    problemData.u1Cont = @(t,x,z) ones(size(x));
    problemData.u2Cont = @(t,x,z) zeros(size(x));
    
    problemData.fhCont = @(t,x) 2 * hVar * x;
    problemData.fuCont = @(t,x,z) 2 * problemData.gConst * hVar * x;
    
    problemData.DCont = { @(t,x,z) ones(size(x)), @(t,x,z) zeros(size(x)); ...
                          @(t,x,z) zeros(size(x)), @(t,x,z) ones(size(x)) };
    
    problemData.hDCont = problemData.hCont;
    problemData.u1DCont = problemData.u1Cont;
    problemData.u2DCont = problemData.u2Cont;
    problemData.q1DCont = @(t,x,z) zeros(size(x));
    problemData.q2DCont = @(t,x,z) zeros(size(x));
    problemData.uhDCont = @(t,x) problemData.hCont(t,x);
    
  case 'quadratic h, linear u'
    domainWidth = 1;
    idLand = -1; idOS = -1; idRiv = -1; idRad = -1;
    
    problemData.gConst = 10;
    xi0Cont = @(x) ones(size(x));
    zBotCont = @(x) zeros(size(x));
    hVar = 0.05;
    
    problemData.hCont = @(t,x) 1 + (x.^2 - domainWidth/2) * hVar;
    problemData.u1Cont = @(t,x,z) x;
    problemData.u2Cont = @(t,x,z) -z;
    
    problemData.fhCont = @(t,x) 2 * hVar * x.^2 + problemData.hCont(t,x);
    problemData.fuCont = @(t,x,z) x + 2 * problemData.gConst * hVar * x;
    
    problemData.DCont = { @(t,x,z) ones(size(x)), @(t,x,z) zeros(size(x)); ...
                          @(t,x,z) zeros(size(x)), @(t,x,z) ones(size(x)) };
    
    problemData.hDCont = problemData.hCont;
    problemData.u1DCont = problemData.u1Cont;
    problemData.u2DCont = problemData.u2Cont;
    problemData.q1DCont = @(t,x,z) -ones(size(x));
    problemData.q2DCont = @(t,x,z) zeros(size(x));
    problemData.uhDCont = @(t,x) x .* problemData.hCont(t,x);

  case 'quadratic u'
    domainWidth = 1;
    idLand = -1; idOS = -1; idRiv = -1; idRad = -1;
    
    problemData.gConst = 10;
    xi0Cont = @(x) ones(size(x));
    zBotCont = @(x) zeros(size(x));
    
    problemData.hCont = @(t,x) ones(size(x));
    problemData.u1Cont = @(t,x,z) x - z.^2;
    problemData.u2Cont = @(t,x,z) -z;
    
    problemData.fhCont = @(t,x) ones(size(x));
    problemData.fuCont = @(t,x,z) x + z.^2 + 2;
    
    problemData.DCont = { @(t,x,z) ones(size(x)), @(t,x,z) zeros(size(x)); ...
                          @(t,x,z) zeros(size(x)), @(t,x,z) ones(size(x)) };
    
    problemData.hDCont = problemData.hCont;
    problemData.u1DCont = problemData.u1Cont;
    problemData.u2DCont = problemData.u2Cont;
    problemData.q1DCont = @(t,x,z) -ones(size(x));
    problemData.q2DCont = @(t,x,z) 2 * z;
    problemData.uhDCont = @(t,x) x - 1/3;

  case 'quadratic'
    domainWidth = 1;
    idLand = -1; idOS = -1; idRiv = -1; idRad = -1;
    
    problemData.gConst = 10;
    xi0Cont = @(x) ones(size(x));
    zBotCont = @(x) zeros(size(x));
    hVar = 0.05;
    
    problemData.hCont = @(t,x) 1 + (x.^2 - domainWidth/2) * hVar;
    problemData.u1Cont = @(t,x,z) x - z.^2;
    problemData.u2Cont = @(t,x,z) -z;
    
    problemData.fhCont = @(t,x) 2 * hVar * x.^2 + problemData.hCont(t,x) .* (1 - 2 * hVar * x .* problemData.hCont(t,x));
    problemData.fuCont = @(t,x,z) (1 + 2 * problemData.gConst * hVar) * x + z.^2 + 2;
    
    problemData.DCont = { @(t,x,z) ones(size(x)), @(t,x,z) zeros(size(x)); ...
                          @(t,x,z) zeros(size(x)), @(t,x,z) ones(size(x)) };
    
    problemData.hDCont = problemData.hCont;
    problemData.u1DCont = problemData.u1Cont;
    problemData.u2DCont = problemData.u2Cont;
    problemData.q1DCont = @(t,x,z) -ones(size(x));
    problemData.q2DCont = @(t,x,z) 2 * z;
    problemData.uhDCont = @(t,x) x .* problemData.hCont(t,x) - problemData.hCont(t,x).^3 / 3;

  case 'linear h, quadratic u, linear D'
    domainWidth = 1;
    idLand = -1; idOS = -1; idRiv = -1; idRad = -1;
    
    problemData.gConst = 10;
    xi0Cont = @(x) ones(size(x));
    zBotCont = @(x) zeros(size(x));
    hVar = 0.05;
    
    problemData.hCont = @(t,x) 1 + (x - domainWidth/2) * hVar;
    problemData.u1Cont = @(t,x,z) x - z.^2;
    problemData.u2Cont = @(t,x,z) -z;
    
    problemData.fhCont = @(t,x) hVar * x + problemData.hCont(t,x) .* (1 - hVar * problemData.hCont(t,x));
    problemData.fuCont = @(t,x,z) x + z.^2 + 4 * z + 1;
    
    problemData.DCont = { @(t,x,z) x + 1,         @(t,x,z) ones(size(x)); ...
                          @(t,x,z) ones(size(x)), @(t,x,z) z + 1          };
    
    problemData.hDCont = problemData.hCont;
    problemData.u1DCont = problemData.u1Cont;
    problemData.u2DCont = problemData.u2Cont;
    problemData.q1DCont = @(t,x,z) -x + 2 * z - 1;
    problemData.q2DCont = @(t,x,z) 2 * z .* (z + 1) - 1;
    problemData.uhDCont = @(t,x) x .* problemData.hCont(t,x) - problemData.hCont(t,x).^3 / 3;

  case 'convergence'
    if license('checkout', 'Symbolic_Toolbox')
      syms x z t

      domainWidth = 100;
      idLand = -1; idOS = -1; idRiv = 2; idRad = 4;

      gSym = sym('10');
      zBotSym(x) = sym('0');

      deltaSym = sym('0.1');
      rhoSym = sym('0.2');

      hSym(t,x) = deltaSym * sin(rhoSym * (t + x)) + 2;
      u1Sym(t,x,z) = sqrt(deltaSym) * z * sin(rhoSym * (t + x));
      u2Sym(t,x,z) = -0.5 * sqrt(deltaSym) * rhoSym * z^2 * cos(rhoSym * (t + x));
      DSym = cellfun(@(c) symfun(c, [t x z]), {0.001, 0; 0, 0.001}, 'UniformOutput', false);

      [problemData, xi0Cont, zBotCont] = analyticalData(problemData, hSym, u1Sym, u2Sym, gSym, zBotSym, DSym, domainWidth);
    else
      error('Symbolic Toolbox required to derive problem formulation!')
    end % if
    
  case 'convergence2'
    if license('checkout', 'Symbolic_Toolbox')
      syms x z t

      domainWidth = 100;
      idLand = -1; idOS = -1; idRiv = -1; idRad = -1;

      gSym = sym('9.81');
      zBotSym(x) = -2 + 0.005 * x;

      epsSym = sym('0.01');
      deltaSym = sym('0.1');
      omegaSym = sym('0.01');

      hSym(t,x) = epsSym * sin(omegaSym * x + t) - zBotSym(x);
      u1Sym(t,x,z) = deltaSym * (z - zBotSym(x)) * sin(omegaSym * x + t);
      u2Sym(t,x,z) = deltaSym * (z - zBotSym(x)) * ( diff(zBotSym, x) * sin(omegaSym * x + t) - omegaSym/2 * (z - zBotSym(x)) * cos(omegaSym * x + t) );
      DSym = cellfun(@(c) symfun(c, [t x z]), {0.1, 0; 0, 0.1}, 'UniformOutput', false);

      [problemData, xi0Cont, zBotCont] = analyticalData(problemData, hSym, u1Sym, u2Sym, gSym, zBotSym, DSym, domainWidth);
    else
      error('Symbolic Toolbox required to derive problem formulation!')
    end % if
end % switch
end % function


function [problemData, xi0Cont, zBotCont] = analyticalData(problemData, hSym, u1Sym, u2Sym, gConst, zBotSym, DSym, domainWidth)
syms x z t

%% Partial derivatives of solution
dxU1Sym = diff(u1Sym, x);
dzU1Sym = diff(u1Sym, z);
dzU2Sym = diff(u2Sym, z);

%% Check continuity
%assert(isequal(dxU1Sym + dzU2Sym, symfun(0, [t x z])), 'u1 and u2 do not fulfill continuity equation')

%% Depth integrated velocity
depthIntU1Sym = int(u1Sym, z, zBotSym(x), zBotSym(x) + hSym(t,x));

%% Compute boundary conditions
q1DSym = DSym{1,1} * dxU1Sym + DSym{1,2} * dzU1Sym;
q2DSym = DSym{2,1} * dxU1Sym + DSym{2,2} * dzU1Sym;

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

problemData.q1DCont = matlabFunction(q1DSym, 'Vars', [t x z]);
problemData.q2DCont = matlabFunction(q2DSym, 'Vars', [t x z]);
problemData.uhDCont = matlabFunction(depthIntU1Sym, 'Vars', [t x]);

%% Determine constants
problemData.gConst = double(gConst);
zBot = matlabFunction(zBotSym, 'Vars', x);
zBotCont = @(x) zBot(x) .* ones(size(x));
xi0Sym = int(hSym(0, x) + zBotSym(x), x, 0, domainWidth) / domainWidth;
xi0Cont = @(x) double(xi0Sym) * ones(size(x));
end % function
