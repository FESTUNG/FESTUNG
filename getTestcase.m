% Define testcases (initial and boundary data, coefficients, etc.) for swe_2dv

%===============================================================================
%> @file
%>
%> @brief Define testcases (initial and boundary data, coefficients, etc.) 
%>        for @link swe_2dv @endlink
%===============================================================================
%>
%> @brief Define testcases (initial and boundary data, coefficients, etc.) 
%>        for @link swe_2dv @endlink
%>
%> This routine is called by @link swe_2dv/configureProblem.m @endlink
%>
%> Currently, the following testcases are available:
%>
%> - showcase
%> - coupled_constXi
%> - coupled_stationary
%> - coupled_transient
%> - dambreak
%> - convergence
%>
%> Please refer to the function definitions in the code for their definitions.
%>
%> A detailed description can be found in @ref RRAFK2018.
%>
%> @param  problemData  A struct with problem parameters as provided by 
%>                      swe_2dv/configureProblem.m . @f$[\text{struct}]@f$
%> @param  problemName  The name of the testcase to be returned.
%>
%> @retval problemData  A struct with all necessary parameters and definitions
%>                      for the problem description. @f$[\text{struct}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2018 Balthasar Reuter, Florian Frank, Vadym Aizinger
%>
%> @author Balthasar Reuter, 2018
%> 
%> @par License
%> @parblock
%> This program is free software: you can redistribute it and/or modify
%> it under the terms of the GNU General Public License as published by
%> the Free Software Foundation, either version 3 of the License, or
%> (at your option) any later version.
%>
%> This program is distributed in the hope that it will be useful,
%> but WITHOUT ANY WARRANTY; without even the implied warranty of
%> MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%> GNU General Public License for more details.
%>
%> You should have received a copy of the GNU General Public License
%> along with this program.  If not, see <http://www.gnu.org/licenses/>.
%> @endparblock
%
function problemData = getTestcase(problemData, problemName)

isHotstart = false;

switch problemName
  case 'showcase'
    isAnalytical = false;
    isHotstart = false;
    hotstartFile = ['swe_2dv' filesep 'showcase_p1_50x10.mat'];
    
    domainWidth = 100;
    idLand = -1; idRad = -1; idRiv = 4; idOS = 2; idBdrRiem = [2,4];
    
    gConst = 10;
    CfConst = 4.e-3;
    rho = 8.e-2;
    
    zBotCont = @(x) (cos((x-50)/30 * pi) + 1) .* (20 <= x & x <= 80);
    xi0Cont = @(x) 4.75 * ones(size(x));
    h0Cont = @(x) 5 - zBotCont(x);
    u10Cont = @(x,z) zeros(size(x));
    
    fhCont = @(t,x) zeros(size(x));
    fuCont = @(t,x,z) zeros(size(x));
    
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    
    u1DCont = @(t,x,z) log(1 + (z - zBotCont(x)) ./ (5 - zBotCont(x)) * (exp(1) - 1));
    hDCont = @(t,x) 5 * ones(size(x));
    u2DCont = @(t,x,z) zeros(size(x));
    q1DCont = @(t,x,z) zeros(size(x));
    q2DCont = @(t,x,z) zeros(size(x));
    
  
  case {'coupled_constXi', 'coupled_stationary', 'coupled_transient'}
    isAnalytical = true;
    
    domainWidth = 100;
    idBdrU = [2, 4]; idBdrH = [2, 4]; idBdrQ = [2, 4]; idBdrRiem = [2, 4];
  
    gConst = 10;
    CfConst = 0;
    
    aConst = 0;  % zBot level
    bConst = 0.005;  % zBot slope
    cConst = 1;  % constant for u1
    dConst = 0.05; % diffusion coefficnet
    alphaConst = 0.1; % amplitude for x2-dependent part of u2
    betaConst = 0.3; % frequency for x2-dependent part of hydraulic head
    gammaConst = 1; % amplitude for x1-dependent part of u1
    deltaConst = 0.07; % frequency for x1-dependent part of u1
    etaConst = 0.003;
    thetaConst = 0.4;
    rhoConst = 0.08;
    tauConst = 0.08;
    kConst = 0.01; % diffusion coefficient for subsurface problem
    kappaConst = 1; % amplitude for x1-dependent part of hydraulic head
    lambdaConst = 0.07; % frequency for x1-dependent part of hydraulic head
    muConst = 0.07;
        
    xi0Cont = @(x) 5 * ones(size(x));
    zBotCont = @(x) aConst + bConst * x;
    
    switch problemName
      case 'coupled_constXi'
        xiCont = @(t,x) xi0Cont(x);
        hCont = @(t,x) xiCont(t,x) - zBotCont(x);

        dxXiCont = @(t,x) zeros(size(x));
        dtHCont = @(t,x) zeros(size(x));
      case 'coupled_stationary'
        xiCont = @(t,x) 5 + etaConst * sin(rhoConst * x);
        hCont = @(t,x) xiCont(t,x) - zBotCont(x);

        dxXiCont = @(t,x) etaConst * rhoConst * cos(rhoConst * x);
        dtHCont = @(t,x) zeros(size(x));
      case 'coupled_transient'
        xiCont = @(t,x) 5 + etaConst * sin(rhoConst * x + tauConst * t);
        hCont = @(t,x) xiCont(t,x) - zBotCont(x);

        dxXiCont = @(t,x) etaConst * rhoConst * cos(rhoConst * x + tauConst * t);
        dtHCont = @(t,x) etaConst * tauConst * cos(rhoConst * x + tauConst * t);
    end % switch
    
    switch problemName
      case {'coupled_constXi', 'coupled_stationary'}
        yCont = @(t,x) gammaConst * sin(deltaConst * x) + cConst;
        dtYCont = @(t,x) zeros(size(x));
        dxYCont = @(t,x) gammaConst * deltaConst * cos(deltaConst * x);
        dxdxYCont = @(t,x) -gammaConst * deltaConst * deltaConst * sin(deltaConst * x);

        omegaCont = @(t,x) kappaConst * cos(lambdaConst * x);
      case {'coupled_transient'}
        yCont = @(t,x) gammaConst * sin(deltaConst * x + thetaConst * t) + cConst;
        dtYCont = @(t,x) gammaConst * thetaConst * cos(deltaConst * x + thetaConst * t);
        dxYCont = @(t,x) gammaConst * deltaConst * cos(deltaConst * x + thetaConst * t);
        dxdxYCont = @(t,x) -gammaConst * deltaConst * deltaConst * sin(deltaConst * x + thetaConst * t);

        omegaCont = @(t,x) kappaConst * cos(lambdaConst * x + muConst * t);
    end % switch
    
    q1ZbCont = @(t,x) -dxXiCont(t,x) + betaConst * bConst * cos(betaConst * zBotCont(x)) .* omegaCont(t,x);
    q2ZbCont = @(t,x) -betaConst * cos(betaConst * zBotCont(x)) .* omegaCont(t,x);
    vCont = @(t,x,z) -dxYCont(t,x) .* ( sin(alphaConst * z) / alphaConst - cos(alphaConst * zBotCont(x)) .* z ) ...
      - alphaConst * bConst * yCont(t,x) .* sin(alphaConst * zBotCont(x)) .* z;
    epsCont = @(t,x) kConst * (-bConst * q1ZbCont(t,x) + q2ZbCont(t,x)) - vCont(t,x,zBotCont(x));
    
    u1Cont = @(t,x,z) yCont(t,x) .* ( cos(alphaConst * z) - cos(alphaConst * zBotCont(x)) );
    u2Cont = @(t,x,z) vCont(t,x,z) + epsCont(t,x);
    
    dtU1Cont = @(t,x,z) dtYCont(t,x) .* ( cos(alphaConst * z) - cos(alphaConst * zBotCont(x)) );
    dxU1Cont = @(t,x,z) dxYCont(t,x) .* ( cos(alphaConst * z) - cos(alphaConst * zBotCont(x)) ) ...
      + alphaConst * bConst * yCont(t,x) .* sin(alphaConst * zBotCont(x));
    dzU1Cont = @(t,x,z) -alphaConst * yCont(t,x) .* sin(alphaConst * z);
    dzU2Cont = @(t,x,z) -dxYCont(t,x) .* ( cos(alphaConst * z) - cos(alphaConst * zBotCont(x)) ) ...
      - alphaConst * bConst * yCont(t,x) .* sin(alphaConst * zBotCont(x));
    
    dxdxU1Cont = @(t,x,z) dxdxYCont(t,x) .* ( cos(alphaConst * z) - cos(alphaConst * zBotCont(x)) ) ...
      + 2 * alphaConst * bConst * dxYCont(t,x) .* sin(alphaConst * zBotCont(x)) ...
      + alphaConst * alphaConst * bConst * bConst * yCont(t,x) .* cos(alphaConst * zBotCont(x));
    dxdzU1Cont = @(t,x,z) -alphaConst * dxYCont(t,x) .* sin(alphaConst * z);
    dzdzU1Cont = @(t,x,z) -alphaConst * alphaConst * yCont(t,x) .* cos(alphaConst * z);
    
    dxU1zIntCont = @(t,x) dxYCont(t,x) .* ( 1/alphaConst * ( sin(alphaConst * xiCont(t,x)) - sin(alphaConst * zBotCont(x)) ) ...
          - cos(alphaConst * zBotCont(x)) .* hCont(t,x) ) ...
        + yCont(t,x) .* ( alphaConst * bConst * sin(alphaConst * zBotCont(x)) .* hCont(t,x) ...
          + ( cos(alphaConst * xiCont(t,x)) - cos(alphaConst * zBotCont(x)) ) .* dxXiCont(t,x) );
    
    DCont = @(t,x,z) dConst * ones(size(x));
    dxzDCont = @(t,x,z) zeros(size(x));
    
  case 'dambreak'
    isAnalytical = false;
    
    domainWidth = 100;
    idLand = -1; idRad = [2,4]; idRiv = -1; idOS = -1; idBdrRiem = -1;
    
    gConst = 10;
    rho = 0.001;
    
    zBotCont = @(x) zeros(size(x));
    xi0Cont = @(x) 2 * ones(size(x));
    h0Cont = @(x) .1 + 4 * (x <= 20);
    u10Cont = @(x,z) zeros(size(x));
    
    fhCont = @(t,x) zeros(size(x));
    fuCont = @(t,x,z) zeros(size(x));
    
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    
    hDCont = @(t,x) zeros(size(x));
    u1DCont = @(t,x,z) zeros(size(x));
    u2DCont = @(t,x,z) zeros(size(x));
    q1DCont = @(t,x,z) zeros(size(x));
    q2DCont = @(t,x,z) zeros(size(x));
  
  case 'convergence'
    isAnalytical = true;
    
    domainWidth = 100;
    idBdrU = [2, 4]; idBdrH = [2, 4]; idBdrQ = [2, 4]; idBdrRiem = [2, 4];
    
    gConst = 10;
    xi0Cont = @(x) zeros(size(x));
    dxZb = -0.005;
    omega = 0.1;
    delta = 0.1;
    epsilon = 0.01;
    rho = 0;%0.001;
    CfConst = 0;
    
    xiCont = @(t,x) epsilon * sin(omega * (x+t));
    zBotCont = @(x) -2 + dxZb * x;
    
    hCont = @(t,x) xiCont(t,x) - zBotCont(x);
    u1Cont = @(t,x,z) delta * (z - zBotCont(x)) .* sin(omega * (x+t));
    u2Cont = @(t,x,z) delta * dxZb * z .* sin(omega * (x+t)) - 0.5 * delta * omega * (z - zBotCont(x)).^2 .* cos(omega * (x+t));
    
    dxXiCont = @(t,x) epsilon * omega * cos(omega * (x+t));
    dtHCont = @(t,x) epsilon * omega * cos(omega * (x+t));
    
    dtU1Cont = @(t,x,z) delta * omega * (z - zBotCont(x)) .* cos(omega * (x+t));
    dxU1Cont = @(t,x,z) -delta * dxZb * sin(omega * (x+t)) + delta * omega * (z - zBotCont(x)) .* cos(omega * (x+t));
    dzU1Cont = @(t,x,z) delta * sin(omega * (x+t));
    dzU2Cont = @(t,x,z) delta * dxZb * sin(omega * (x+t)) - delta * omega * (z - zBotCont(x)) .* cos(omega * (x+t));
    
    dxdxU1Cont = @(t,x,z) -delta * omega * ( 2 * dxZb * cos(omega * (x+t)) + omega * (z - zBotCont(x)) .* sin(omega * (x+t)) );
    dxdzU1Cont = @(t,x,z) delta * omega * cos(omega * (x+t));
    dzdzU1Cont = @(t,x,z) zeros(size(x));
    
    dxU1zIntCont = @(t,x) 0.5 * delta * omega * cos(omega * (x+t)) .* (xiCont(t,x) - zBotCont(x)).^2 + ...
      delta * sin(omega * (x+t)) .* (xiCont(t,x) - zBotCont(x)) .* (dxXiCont(t,x) - dxZb);
    
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
      @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    dxzDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
      @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
        
  otherwise
    error('ERROR: unknown testcase')
end % switch

problemData = setdefault(problemData, 'isHotstart', isHotstart);
if problemData.isHotstart
  problemData = setdefault(problemData, 'hotstartFile', hotstartFile);
end % if

problemData = setdefault(problemData, 'domainWidth', domainWidth);
problemData = setdefault(problemData, 'gConst', gConst);
problemData = setdefault(problemData, 'CfConst', CfConst);
problemData = setdefault(problemData, 'xi0Cont', xi0Cont);
problemData = setdefault(problemData, 'zBotCont', zBotCont);

% Convert idLand, idRiver, idOS to idBdrU, idBdrQ, idBdrH
if ~isequal([exist('idBdrU', 'var'), exist('idBdrH', 'var'), exist('idBdrQ', 'var')], [1,1,1])
  assert(isequal([exist('idLand', 'var'), exist('idRiv', 'var'), exist('idOS', 'var'), exist('idRad', 'var')], [1,1,1,1]), 'No boundary conditions given.')
  assert(~any(intersect(idLand(:), idRad(:)) > 0), 'Land and Radiation on same boundary specified.')
  assert(~any(intersect(idLand(:), idRiv(:)) > 0), 'Land and River on same boundary specified.')
  assert(~any(intersect(idLand(:), idOS(:)) > 0), 'Land and Open Sea on same boundary specified.')
  assert(~any(intersect(idRad(:), idRiv(:)) > 0), 'Radiation and River on same boundary specified.')
  assert(~any(intersect(idRad(:), idOS(:)) > 0), 'Radiation and Open Sea on same boundary specified.')
  assert(~any(intersect(idRiv(:), idOS(:)) > 0), 'River and Open Sea on same boundary specified.')
  fprintf('Converting physical to numerical boundary conditions.\n')
  idBdrU = sort([idLand, idRiv]);
  idBdrH = sort([idRiv, idOS]);
  idBdrQ = sort([idRad, idOS]);
end % if

problemData = setdefault(problemData, 'idBdrRiem', idBdrRiem);
problemData = setdefault(problemData, 'idBdrH', idBdrH);
problemData = setdefault(problemData, 'idBdrU', idBdrU);
problemData = setdefault(problemData, 'idBdrQ', idBdrQ);

if isAnalytical
  t0 = problemData.t0;
  
  problemData = setdefault(problemData, 'hCont', hCont);
  problemData = setdefault(problemData, 'u1Cont', u1Cont);
  problemData = setdefault(problemData, 'u2Cont', u2Cont);
  problemData = setdefault(problemData, 'DCont', DCont);

  fhCont = @(t,x) dtHCont(t,x) + dxU1zIntCont(t,x);
  
  if iscell(DCont) && ~isvector(DCont)
    fuCont = @(t,x,z) dtU1Cont(t,x,z) + 2 * u1Cont(t,x,z) .* dxU1Cont(t,x,z) + ...
      dzU1Cont(t,x,z) .* u2Cont(t,x,z) + u1Cont(t,x,z) .* dzU2Cont(t,x,z) - ( ...
      DCont{1,1}(t,x,z) .* dxdxU1Cont(t,x,z) + dxzDCont{1,1}(t,x,z) .* dxU1Cont(t,x,z) + ...
      DCont{1,2}(t,x,z) .* dxdzU1Cont(t,x,z) + dxzDCont{1,2}(t,x,z) .* dzU1Cont(t,x,z) + ...
      DCont{2,1}(t,x,z) .* dxdzU1Cont(t,x,z) + dxzDCont{2,1}(t,x,z) .* dxU1Cont(t,x,z) + ...
      DCont{2,2}(t,x,z) .* dzdzU1Cont(t,x,z) + dxzDCont{2,2}(t,x,z) .* dzU1Cont(t,x,z) ) + ...
      gConst * dxXiCont(t,x);
    q1DCont = @(t,x,z) -DCont{1,1}(t,x,z) .* dxU1Cont(t,x,z) - DCont{1,2}(t,x,z) .* dzU1Cont(t,x,z);
    q2DCont = @(t,x,z) -DCont{2,1}(t,x,z) .* dxU1Cont(t,x,z) - DCont{2,2}(t,x,z) .* dzU1Cont(t,x,z);
  elseif iscell(DCont) && isvector(DCont)
    fuCont = @(t,x,z) dtU1Cont(t,x,z) + 2 * u1Cont(t,x,z) .* dxU1Cont(t,x,z) + ...
      dzU1Cont(t,x,z) .* u2Cont(t,x,z) + u1Cont(t,x,z) .* dzU2Cont(t,x,z) - ( ...
      DCont{1}(t,x,z) .* dxdxU1Cont(t,x,z) + dxzDCont{1}(t,x,z) .* dxU1Cont(t,x,z) + ...
      DCont{2}(t,x,z) .* dzdzU1Cont(t,x,z) + dxzDCont{2}(t,x,z) .* dzU1Cont(t,x,z) ) + ...
      gConst * dxXiCont(t,x);
    q1DCont = @(t,x,z) -DCont{1}(t,x,z) .* dxU1Cont(t,x,z);
    q2DCont = @(t,x,z) -DCont{2}(t,x,z) .* dzU1Cont(t,x,z);
  else
    fuCont = @(t,x,z) dtU1Cont(t,x,z) + 2 * u1Cont(t,x,z) .* dxU1Cont(t,x,z) + ...
      dzU1Cont(t,x,z) .* u2Cont(t,x,z) + u1Cont(t,x,z) .* dzU2Cont(t,x,z) - ( ...
      DCont(t,x,z) .* dxdxU1Cont(t,x,z) + dxzDCont(t,x,z) .* dxU1Cont(t,x,z) + ...
      DCont(t,x,z) .* dzdzU1Cont(t,x,z) + dxzDCont(t,x,z) .* dzU1Cont(t,x,z) ) + ...
      gConst * dxXiCont(t,x);
    q1DCont = @(t,x,z) -DCont(t,x,z) .* dxU1Cont(t,x,z);
    q2DCont = @(t,x,z) -DCont(t,x,z) .* dzU1Cont(t,x,z);
  end % if

  problemData = setdefault(problemData, 'fhCont', fhCont);
  problemData = setdefault(problemData, 'fuCont', fuCont);

  problemData = setdefault(problemData, 'h0Cont', @(x1) hCont(t0, x1));
  problemData = setdefault(problemData, 'u10Cont', @(x1,x2) u1Cont(t0, x1, x2));

  problemData = setdefault(problemData, 'hDCont', hCont);
  problemData = setdefault(problemData, 'u1DCont', u1Cont);
  problemData = setdefault(problemData, 'u2DCont', u2Cont);
  problemData = setdefault(problemData, 'q1DCont', q1DCont);
  problemData = setdefault(problemData, 'q2DCont', q2DCont);
else
  problemData = setdefault(problemData, 'fhCont', fhCont);
  problemData = setdefault(problemData, 'fuCont', fuCont);
  problemData = setdefault(problemData, 'DCont', DCont);
  
  problemData = setdefault(problemData, 'h0Cont', h0Cont);
  problemData = setdefault(problemData, 'u10Cont', u10Cont);

  problemData = setdefault(problemData, 'hDCont', hDCont);
  problemData = setdefault(problemData, 'u1DCont', u1DCont);
  problemData = setdefault(problemData, 'u2DCont', u2DCont);
  problemData = setdefault(problemData, 'q1DCont', q1DCont);
  problemData = setdefault(problemData, 'q2DCont', q2DCont);
end % if
end % function
