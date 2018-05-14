function problemData = getTestcase(problemData, problemName)

isAnalytical = true;
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
    
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
      @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    
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
    
    aConst = 0;
    bConst = 0.005;
    cConst = 1;
    dConst = 0.05;
    alphaConst = 0.1;
    betaConst = 0.3;
    gammaConst = 1;
    deltaConst = 0.07;
    etaConst = 0.003;
    thetaConst = 0.4;
    rhoConst = 0.08;
    tauConst = 0.08;
    kConst = 0.01;
    kappaConst = 1;
    lambdaConst = 0.07;
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
    
    DCont = cellfun(@(coef) @(t,x,z) coef * ones(size(x)), {dConst, 0; 0, dConst}, 'UniformOutput', false);
    dxzDCont = cellfun(@(coef) @(t,x,z) coef * ones(size(x)), {0, 0; 0, 0}, 'UniformOutput', false);
    
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
    
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
      @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    
    hDCont = @(t,x) zeros(size(x));
    u1DCont = @(t,x,z) zeros(size(x));
    u2DCont = @(t,x,z) zeros(size(x));
    q1DCont = @(t,x,z) zeros(size(x));
    q2DCont = @(t,x,z) zeros(size(x));
  
  case 'convergence'
    domainWidth = 100;
    idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idBdrRiem = [2,4];
    %     idLand = -1; idOS = -1; idRiv = -1; idRad = -1; idRiem = -1;
    %     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = [2,4]; idRiem = -1;
    
    gConst = 10;
    xi0Cont = @(x) zeros(size(x));
    dxZb = -0.005;
    omega = 0.1;
    delta = 0.1;
    epsilon = 0.01;
    rho = 0;%0.001;
    
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
    
  case 'coupling'
    domainWidth = 100;
    idBdrU = [2, 4]; idBdrH = [2, 4]; idBdrQ = [2, 4]; idBdrRiem = [2, 4];
    
    gConst = 10;
    xi0Cont = @(x) 5 * ones(size(x));
    dxZb = 0;
    a = 0.05;
    b = 0.1;
    c = 0.1;
    d = 0.01;
    k = 1;
    
    xiCont = @(t,x) 5 + a * cos(b*x + c*t);
    zBotCont = @(x) 0 + dxZb * x;
    
    hCont = @(t,x) xiCont(t,x) - zBotCont(x);
    u1Cont = @(t,x,z) a * b * k * sin(b*x + c*t) .* exp(-b*z);
    u2Cont = @(t,x,z) a * b * k * cos(b*x + c*t) .* exp(-b*z);
    
    dxXiCont = @(t,x) -a * b * sin(b*x + c*t);
    dtHCont = @(t,x) -a * c * sin(b*x + c*t);
    
    dtU1Cont = @(t,x,z) a * b * c * k * cos(b*x + c*t) .* exp(-b*z);
    dxU1Cont = @(t,x,z) a * b^2 * k * cos(b*x + c*t) .* exp(-b*z);
    dzU1Cont = @(t,x,z) -a * b^2 * k * sin(b*x + c*t) .* exp(-b*z);
    dzU2Cont = @(t,x,z) -a * b^2 * k * cos(b*x + c*t) .* exp(-b*z);
    
    dxdxU1Cont = @(t,x,z) -a * b^3 * k * sin(b*x + c*t) .* exp(-b*z);
    dxdzU1Cont = @(t,x,z) -a * b^3 * k * cos(b*x + c*t) .* exp(-b*z);
    dzdzU1Cont = @(t,x,z) a * b^3 * k * sin(b*x + c*t) .* exp(-b*z);
    
    dxU1zIntCont = @(t,x) a * b * k * ( ...
      sin(b*x + c*t) .* (exp(-b * xiCont(t,x)) .* dxXiCont(t,x) - exp(-b * zBotCont(x)) .* dxZb) - ...
      cos(b*x + c*t) .* (exp(-b * xiCont(t,x)) - exp(-b * zBotCont(x))) );
    
    DCont = cellfun(@(coef) @(t,x,z) coef * ones(size(x)), {d, 0; 0, d}, 'UniformOutput', false);
    dxzDCont = cellfun(@(coef) @(t,x,z) coef * ones(size(x)), {0, 0; 0, 0}, 'UniformOutput', false);
    
  case 'coupling2'
    domainWidth = 100;
    idBdrU = [2, 4]; idBdrH = [2, 4]; idBdrQ = [2, 4]; idBdrRiem = [2, 4];
    
    gConst = 10;
    xi0Cont = @(x) 5 * ones(size(x));
    dxZb = 0.005;
    a = 0.05;
    b = 0.1;
    c = 0.1;
    d = 0.01;
    k = 1;
    
    zBotCont = @(x) 0 + dxZb * x;
    xiCont = @(t,x) 5 + a * cos(b*x + c*t) .* exp(-b * dxZb * x) + zBotCont(x);
    
    hCont = @(t,x) xiCont(t,x) - zBotCont(x);
    u1Cont = @(t,x,z) a * b * k * sin(b*x + c*t) .* exp(-b*z);
    u2Cont = @(t,x,z) a * b * k * cos(b*x + c*t) .* exp(-b*z);
    
    dxXiCont = @(t,x) -a * b * sin(b*x + c*t) .* exp(-b * dxZb * x) - a * b * dxZb * cos(b*x + c*t) .* exp(-b * dxZb * x) + dxZb;
    dtHCont = @(t,x) -a * c * sin(b*x + c*t) .* exp(-b * dxZb * x);
    
    dtU1Cont = @(t,x,z) a * b * c * k * cos(b*x + c*t) .* exp(-b*z);
    dxU1Cont = @(t,x,z) a * b^2 * k * cos(b*x + c*t) .* exp(-b*z);
    dzU1Cont = @(t,x,z) -a * b^2 * k * sin(b*x + c*t) .* exp(-b*z);
    dzU2Cont = @(t,x,z) -a * b^2 * k * cos(b*x + c*t) .* exp(-b*z);
    
    dxdxU1Cont = @(t,x,z) -a * b^3 * k * sin(b*x + c*t) .* exp(-b*z);
    dxdzU1Cont = @(t,x,z) -a * b^3 * k * cos(b*x + c*t) .* exp(-b*z);
    dzdzU1Cont = @(t,x,z) a * b^3 * k * sin(b*x + c*t) .* exp(-b*z);
    
    dxU1zIntCont = @(t,x) a * b * k * ( ...
      sin(b*x + c*t) .* (exp(-b * xiCont(t,x)) .* dxXiCont(t,x) - exp(-b * zBotCont(x)) .* dxZb) - ...
      cos(b*x + c*t) .* (exp(-b * xiCont(t,x)) - exp(-b * zBotCont(x))) );
    
    DCont = cellfun(@(coef) @(t,x,z) coef * ones(size(x)), {d, 0; 0, d}, 'UniformOutput', false);
    dxzDCont = cellfun(@(coef) @(t,x,z) coef * ones(size(x)), {0, 0; 0, 0}, 'UniformOutput', false);
    
  case 'coupling_linear'
    domainWidth = 100;
    idBdrU = [2, 4]; idBdrH = [2, 4]; idBdrQ = [2, 4]; idBdrRiem = [2, 4];
    
    gConst = 10;
    xi0Cont = @(x) 5 * ones(size(x));
    dxZb = 0;
    a = 0.01;
    b = 0.1;
    c = 0.01;
    d = 0.001;
    
    xiCont = @(t,x) 5 + a * x;
    zBotCont = @(x) 0 + dxZb * x;
    
    hCont = @(t,x) xiCont(t,x) - zBotCont(x);
    u1Cont = @(t,x,z) -a * x - b * (zBotCont(x) - z);
    u2Cont = @(t,x,z) -a * (zBotCont(x) - z) + b * dxZb * z - c * x;
    
    dxXiCont = @(t,x) a * ones(size(x));
    dtHCont = @(t,x) zeros(size(x));
    
    dtU1Cont = @(t,x,z) zeros(size(x));
    dxU1Cont = @(t,x,z) -(a + b * dxZb) * ones(size(x));
    dzU1Cont = @(t,x,z) b * ones(size(x));
    dzU2Cont = @(t,x,z) (a + b * dxZb) * ones(size(x));
    
    dxdxU1Cont = @(t,x,z) zeros(size(x));
    dxdzU1Cont = @(t,x,z) zeros(size(x));
    dzdzU1Cont = @(t,x,z) zeros(size(x));
    
    dxU1zIntCont = @(t,x) -a * hCont(t,x) - a * x .* (dxXiCont(t,x) - dxZb) - ...
      b * (dxZb - dxXiCont(t,x) / 2) .* xiCont(t,x) - ...
      b * (zBotCont(x) - xiCont(t,x) / 2) .* dxXiCont(t,x) + b * dxZb * zBotCont(x);
    
    DCont = cellfun(@(c) @(t,x,z) c * ones(size(x)), {d, 0; 0, d}, 'UniformOutput', false);
    dxzDCont = cellfun(@(c) @(t,x,z) c * ones(size(x)), {0, 0; 0, 0}, 'UniformOutput', false);
    
  case 'coupling_rupp'
    domainWidth = 100;
    idBdrU = [2, 4]; idBdrH = [2, 4]; idBdrQ = [2, 4]; idBdrRiem = [2, 4];
    
    gConst = 10;
    xi0Cont = @(x) 5 * ones(size(x));
    dxZb = 0;
    a = 0.01;
    b = 0.3;
    c = 0.5;
    d = 2;
    f = 1;
    k = 1;
    
    xiCont = @(t,x) 5 + a * cos(b*x+c*t);
    zBotCont = @(x) 0 + dxZb * x;
    
    hCont = @(t,x) xiCont(t,x) - zBotCont(x);
    u1Cont = @(t,x,z) -k * a * d^2 / b * sin(b*x+c*t) .* cos(d*z) + f * z;
    u2Cont = @(t,x,z) k * a * d * cos(b*x+c*t) .* sin(d*z) + 1;
    
    dxXiCont = @(t,x) -a * b * sin(b*x+c*t);
    dtHCont = @(t,x) -a * c * sin(b*x+c*t);
    
    dtU1Cont = @(t,x,z) -k * a * c * d^2 / b * cos(b*x+c*t) .* cos(d*z);
    dxU1Cont = @(t,x,z) -k * a * d^2 * cos(b*x+c*t) .* cos(d*z);
    dzU1Cont = @(t,x,z) k * a * d^3 / b * sin(b*x+c*t) .* sin(d*z) + f;
    dzU2Cont = @(t,x,z) k * a * d^2 * cos(b*x+c*t) .* cos(d*z);
    
    dxdxU1Cont = @(t,x,z) k * a * b * d^2 * sin(b*x+c*t) .* cos(d*z);
    dxdzU1Cont = @(t,x,z) k * a * d^3 * cos(b*x+c*t) .* sin(d*z);
    dzdzU1Cont = @(t,x,z) k * a * d^4 / b * sin(b*x+c*t) .* cos(d*z);
    
    dxU1zIntCont = @(t,x) -k * a * d * cos(b*x+c*t) .* ( sin(d*xiCont(t,x)) - sin(d*zBotCont(x)) ) - ...
      k * a * d^2 / b * sin(b*x+c*t) .* ( cos(d*xiCont(t,x)) .* dxXiCont(t,x) - cos(d*zBotCont(x)) * dxZb ) + ...
      f * ( xiCont(t,x) .* dxXiCont(t,x) - zBotCont(x) * dxZb );
    
    DCont = cellfun(@(c) @(t,x,z) c * ones(size(x)), {0.01, 0; 0, 0.01}, 'UniformOutput', false);
    dxzDCont = cellfun(@(c) @(t,x,z) c * ones(size(x)), {0, 0; 0, 0}, 'UniformOutput', false);
    
  case 'utbest'
    domainWidth = 100;
    idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idBdrRiem = [2,4];
    %     idLand = -1; idOS = -1; idRiv = -1; idRad = -1;
    
    gConst = 10;
    xi0Cont = @(x) zeros(size(x));
    dxZb = -0.005;
    omega = 0.15;
    d = 0.125;
    e = 0.125;
    D = 0.125;
    
    xiCont = @(t,x) e * sin(omega * (x+t));
    zBotCont = @(x) -2 + dxZb * x;
    
    hCont = @(t,x) xiCont(t,x) - zBotCont(x);
    u1Cont = @(t,x,z) d * (z - zBotCont(x)) .* sin(omega * (x+t));
    u2Cont = @(t,x,z) d * dxZb * z .* sin(omega * (x+t)) - 0.5 * d * omega * (z - zBotCont(x)).^2 .* cos(omega * (x+t));
    
    dxXiCont = @(t,x) e * omega * cos(omega * (x+t));
    dtHCont = @(t,x) e * omega * cos(omega * (x+t));
    
    dtU1Cont = @(t,x,z) d * omega * (z - zBotCont(x)) .* cos(omega * (x+t));
    dxU1Cont = @(t,x,z) -d * dxZb * sin(omega * (x+t)) + d * omega * (z - zBotCont(x)) .* cos(omega * (x+t));
    dzU1Cont = @(t,x,z) d * sin(omega * (x+t));
    dzU2Cont = @(t,x,z) d * dxZb * sin(omega * (x+t)) - d * omega * (z - zBotCont(x)) .* cos(omega * (x+t));
    
    dxdxU1Cont = @(t,x,z) -d * omega * ( 2 * dxZb * cos(omega * (x+t)) + omega * (z - zBotCont(x)) .* sin(omega * (x+t)) );
    dxdzU1Cont = @(t,x,z) d * omega * cos(omega * (x+t));
    dzdzU1Cont = @(t,x,z) zeros(size(x));
    
    dxU1zIntCont = @(t,x) 0.5 * d * omega * cos(omega * (x+t)) .* hCont(t,x).^2 + ...
      d * sin(omega * (x+t)) .* hCont(t,x) .* (dxXiCont(t,x) - dxZb);
    
    DCont = cellfun(@(c) @(t,x,z) c * ones(size(x)), {D, 0; 0, D}, 'UniformOutput', false);
    dxzDCont = cellfun(@(c) @(t,x,z) c * ones(size(x)), {0, 0; 0, 0}, 'UniformOutput', false);
    
  case 'utbest_sinus'
    domainWidth = 100;
    idBdrU = [2, 4]; idBdrH = [2, 4]; idBdrQ = [2, 4]; idBdrRiem = [2, 4];
    
    gConst = 10;
    xi0Cont = @(x) zeros(size(x));
    dxZb = 0.005;
    omega = 0.1;
    d = 0.3;
    e = 0.1;
    t_coef = 1.;
    D = 0.1;
    
    xiCont = @(t,x) .25 + e * sin(omega * (x+t_coef*t));
    zBotCont = @(x) -2 + dxZb * x;
    
    hCont = @(t,x) xiCont(t,x) - zBotCont(x);
    u1Cont = @(t,x,z) d * ( sin(omega * (x+t_coef*t)) + sin(omega * (z+t_coef*t)) );
    u2Cont = @(t,x,z) -d * omega * z .* cos(omega * (x+t_coef*t));
    
    dxXiCont = @(t,x) e * omega * cos(omega * (x+t_coef*t));
    dtHCont = @(t,x)  e * omega * t_coef * cos(omega * (x+t_coef*t));
    
    dtU1Cont = @(t,x,z) d * omega * t_coef * ( cos(omega * (x+t_coef*t)) + cos(omega * (z+t_coef*t)) );
    dxU1Cont = @(t,x,z) d * omega * cos(omega * (x+t_coef*t));
    dzU1Cont = @(t,x,z) d * omega * cos(omega * (z+t_coef*t));
    dzU2Cont = @(t,x,z) -d * omega * cos(omega * (x+t_coef*t));
    
    dxdxU1Cont = @(t,x,z) -d * omega^2 * sin(omega * (x+t_coef*t));
    dxdzU1Cont = @(t,x,z) zeros(size(x));
    dzdzU1Cont = @(t,x,z) -d * omega^2 * sin(omega * (z+t_coef*t));
    
    dxU1zIntCont = @(t,x) d * ( omega * cos(omega * (x+t_coef*t)) .* hCont(t,x) + ...
      sin(omega * (x+t_coef*t)) .* (dxXiCont(t,x) - dxZb) + ...
      sin(omega * (xiCont(t,x)+t_coef*t)) .* dxXiCont(t,x) - ...
      sin(omega * (zBotCont(x)+t_coef*t)) .* dxZb );
    
    DCont = cellfun(@(c) @(t,x,z) c * ones(size(x)), {D, 0; 0, D}, 'UniformOutput', false);
    dxzDCont = cellfun(@(c) @(t,x,z) c * ones(size(x)), {0, 0; 0, 0}, 'UniformOutput', false);
    
  case 'constant' % OK
    domainWidth = 100;
    gConst = 10;
    xi0Cont = @(x) zeros(size(x));
    idBdrU = -1; idBdrH = [2, 4]; idBdrQ = -1; idBdrRiem = [2, 4];
    
    dxZb = -0.005;
    dtXi = 0*0.01;
    rho = 0;
    
    xiCont = @(t,x) dtXi * t * ones(size(x));
    zBotCont = @(x) -2 + dxZb * x;
    
    hCont = @(t,x) xiCont(t,x) - zBotCont(x);
    u1Cont = @(t,x,z) 0.1 * ones(size(x));
    u2Cont = @(t,x,z) zeros(size(x));
    
    dxXiCont = @(t,x) zeros(size(x));
    dtHCont = @(t,x) dtXi * ones(size(x));
    
    dtU1Cont = @(t,x,z) zeros(size(x));
    dxU1Cont = @(t,x,z) zeros(size(x));
    dzU1Cont = @(t,x,z) zeros(size(x));
    dzU2Cont = @(t,x,z) zeros(size(x));
    
    dxdxU1Cont = @(t,x,z) zeros(size(x));
    dxdzU1Cont = @(t,x,z) zeros(size(x));
    dzdzU1Cont = @(t,x,z) zeros(size(x));
    
    dxU1zIntCont = @(t,x) 0.1 * (dxXiCont(t,x) - dxZb);
    
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
      @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    dxzDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
      @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
    
  case 'linear_h' % OK
    domainWidth = 100;
    idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idBdrRiem = [2,4];
    %     idLand = -1; idOS = -1; idRiv = -1; idRad = -1; idRiem = -1;
    %     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = -1;
    
    gConst = 10;
    xi0Cont = @(x) zeros(size(x));
    dxZb = -0.005;
    dxXi = 0.005;
    rho = 0;
    
    xiCont = @(t,x) dxXi * x;
    zBotCont = @(x) -2 + dxZb * x;
    
    hCont = @(t,x) xiCont(t,x) - zBotCont(x);
    u1Cont = @(t,x,z) 0.1 * ones(size(x));
    u2Cont = @(t,x,z) zeros(size(x));
    
    dxXiCont = @(t,x) dxXi * ones(size(x));
    dtHCont = @(t,x) zeros(size(x));
    
    dtU1Cont = @(t,x,z) zeros(size(x));
    dxU1Cont = @(t,x,z) zeros(size(x));
    dzU1Cont = @(t,x,z) zeros(size(x));
    dzU2Cont = @(t,x,z) zeros(size(x));
    
    dxdxU1Cont = @(t,x,z) zeros(size(x));
    dxdzU1Cont = @(t,x,z) zeros(size(x));
    dzdzU1Cont = @(t,x,z) zeros(size(x));
    
    dxU1zIntCont = @(t,x) 0.1 * (dxXiCont(t,x) - dxZb);
    
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
      @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    dxzDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
      @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
    
  case 'z-linear_u' % OK
    domainWidth = 100;
    idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idBdrRiem = [2,4];
    %     idLand = -1; idOS = -1; idRiv = -1; idRad = -1; idRiem = -1;
    %     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = -1;
    
    gConst = 10;
    xi0Cont = @(x) zeros(size(x));
    dzU1 = 0.1;
    dxZb = -0.005;
    rho = 0;
    
    xiCont = @(t,x) zeros(size(x));
    zBotCont = @(x) -2 + dxZb * x;
    
    hCont = @(t,x) xiCont(t,x) - zBotCont(x);
    u1Cont = @(t,x,z) dzU1 * z;
    u2Cont = @(t,x,z) zeros(size(x));
    
    dxXiCont = @(t,x) zeros(size(x));
    dtHCont = @(t,x) zeros(size(x));
    
    dtU1Cont = @(t,x,z) zeros(size(x));
    dxU1Cont = @(t,x,z) zeros(size(x));
    dzU1Cont = @(t,x,z) dzU1 * ones(size(x));
    dzU2Cont = @(t,x,z) zeros(size(x));
    
    dxdxU1Cont = @(t,x,z) zeros(size(x));
    dxdzU1Cont = @(t,x,z) zeros(size(x));
    dzdzU1Cont = @(t,x,z) zeros(size(x));
    
    dxU1zIntCont = @(t,x) 0.1 * (xiCont(t,x) .* dxXiCont(t,x) - zBotCont(x) .* dxZb);
    
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
      @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    dxzDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
      @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
    
  case 'x-linear_w' % OK
    domainWidth = 100;
    %     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = [2,4];
    %     idLand = -1; idOS = -1; idRiv = -1; idRad = -1; idRiem = -1;
    idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idBdrRiem = -1;
    
    gConst = 10;
    xi0Cont = @(x) zeros(size(x));
    dxU2 = 0.01;
    dxZb = -0.005;
    rho = 0;
    
    xiCont = @(t,x) zeros(size(x));
    zBotCont = @(x) -2 + dxZb * x;
    
    hCont = @(t,x) xiCont(t,x) - zBotCont(x);
    u1Cont = @(t,x,z) 0.1 * ones(size(x));
    u2Cont = @(t,x,z) dxU2 * x;
    
    dxXiCont = @(t,x) zeros(size(x));
    dtHCont = @(t,x) zeros(size(x));
    
    dtU1Cont = @(t,x,z) zeros(size(x));
    dxU1Cont = @(t,x,z) zeros(size(x));
    dzU1Cont = @(t,x,z) zeros(size(x));
    dzU2Cont = @(t,x,z) zeros(size(x));
    
    dxdxU1Cont = @(t,x,z) zeros(size(x));
    dxdzU1Cont = @(t,x,z) zeros(size(x));
    dzdzU1Cont = @(t,x,z) zeros(size(x));
    
    dxU1zIntCont = @(t,x) 0.1 * (dxXiCont(t,x) - dxZb);
    
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
      @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    dxzDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
      @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
    
  case 'linear_vel' % OK
    domainWidth = 100;
    idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idBdrRiem = [2,4];
    %     idLand = -1; idOS = -1; idRiv = -1; idRad = -1; idRiem = -1;
    %     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = -1;
    
    gConst = 10;
    xi0Cont = @(x) zeros(size(x));
    dxU1 = -0.01;
    rho = 0;
    
    xiCont = @(t,x) zeros(size(x));
    zBotCont = @(x) -2 * ones(size(x));
    
    hCont = @(t,x) xiCont(t,x) - zBotCont(x);
    u1Cont = @(t,x,z) dxU1 * x;
    u2Cont = @(t,x,z) -dxU1 * z;
    
    dxXiCont = @(t,x) zeros(size(x));
    dtHCont = @(t,x) zeros(size(x));
    
    dtU1Cont = @(t,x,z) zeros(size(x));
    dxU1Cont = @(t,x,z) dxU1 * ones(size(x));
    dzU1Cont = @(t,x,z) zeros(size(x));
    dzU2Cont = @(t,x,z) -dxU1 * ones(size(x));
    
    dxdxU1Cont = @(t,x,z) zeros(size(x));
    dxdzU1Cont = @(t,x,z) zeros(size(x));
    dzdzU1Cont = @(t,x,z) zeros(size(x));
    
    dxU1zIntCont = @(t,x) 2 * dxU1 * ones(size(x));
    
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
      @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    dxzDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
      @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
    
  case 'linear' % OK
    domainWidth = 100;
    gConst = 10;
    xi0Cont = @(x) zeros(size(x));
    idBdrU = -1; idBdrH = [2, 4]; idBdrQ = -1; idBdrRiem = -[2, 4];
    
    dtXi = 0.01;
    dxZb = -0.005;
    dxXi = 0.005;
    dxU1 = 0.01;
    dzU2 = -dxU1;
    rho = 0.1;
    
    xiCont = @(t,x) dxXi * x + dtXi * t;
    zBotCont = @(x) -2 + dxZb * x;
    
    hCont = @(t,x) xiCont(t,x) - zBotCont(x);
    u1Cont = @(t,x,z) dxU1 * x;
    u2Cont = @(t,x,z) dzU2 * z;
    
    dxXiCont = @(t,x) dxXi * ones(size(x));
    dtHCont = @(t,x) dtXi * ones(size(x));
    
    dtU1Cont = @(t,x,z) zeros(size(x));
    dxU1Cont = @(t,x,z) dxU1 * ones(size(x));
    dzU1Cont = @(t,x,z) zeros(size(x));
    dzU2Cont = @(t,x,z) dzU2 * ones(size(x));
    
    dxdxU1Cont = @(t,x,z) zeros(size(x));
    dxdzU1Cont = @(t,x,z) zeros(size(x));
    dzdzU1Cont = @(t,x,z) zeros(size(x));
    
    dxU1zIntCont = @(t,x) dxU1 * hCont(t,x) + (dxU1 * x) .* (dxXiCont(t,x) - dxZb);
    
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
      @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    dxzDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
      @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
    
  case 'quadratic_h' % not ok
    domainWidth = 100;
    %     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = [2,4];
    %     idLand = -1; idOS = -1; idRiv = -1; idRad = -1; idRiem = -1;
    idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idBdrRiem = -1;
    
    gConst = 40;
    xi0Cont = @(x) zeros(size(x));
    dtXi = 0;%0.01;
    dxZb = 0;
    dxXi = -0.3;
    dxU1 = 0;
    dzU2 = -dxU1;
    rho = 0;
    
    xiCont = @(t,x) dxXi * (x/domainWidth - 0.5).^2 + dtXi * t;
    zBotCont = @(x) -2 + dxZb * x;
    
    hCont = @(t,x) xiCont(t,x) - zBotCont(x);
    u1Cont = @(t,x,z) dxU1 * x;
    u2Cont = @(t,x,z) dzU2 * z;
    
    dxXiCont = @(t,x) 2 * dxXi/domainWidth * (x/domainWidth - 0.5);
    dtHCont = @(t,x) dtXi * ones(size(x));
    
    dtU1Cont = @(t,x,z) zeros(size(x));
    dxU1Cont = @(t,x,z) dxU1 * ones(size(x));
    dzU1Cont = @(t,x,z) zeros(size(x));
    dzU2Cont = @(t,x,z) dzU2 * ones(size(x));
    
    dxdxU1Cont = @(t,x,z) zeros(size(x));
    dxdzU1Cont = @(t,x,z) zeros(size(x));
    dzdzU1Cont = @(t,x,z) zeros(size(x));
    
    dxU1zIntCont = @(t,x) dxU1 * hCont(t,x) + (dxU1 * x) .* (dxXiCont(t,x) - dxZb);
    
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
      @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    dxzDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
      @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
    
  case 'quadratic_u' % ok for p=2
    domainWidth = 100;
    idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idBdrRiem = [2,4];
    %     idLand = -1; idOS = -1; idRiv = -1; idRad = -1; idRiem = -1;
    %     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = -1;
    
    gConst = 10;
    xi0Cont = @(x) zeros(size(x));
    dtXi = 0;%0.01;
    dxZb = 0;%-0.005;
    dxXi = 0;%-0.3;
    dxU1 = 0.01;
    rho = 0;
    
    xiCont = @(t,x) dxXi * (x/domainWidth - 0.5).^2 + dtXi * t;
    zBotCont = @(x) -2 + dxZb * x;
    
    hCont = @(t,x) xiCont(t,x) - zBotCont(x);
    u1Cont = @(t,x,z) dxU1 * (x/domainWidth - 0.5).^2 + 0.01;
    u2Cont = @(t,x,z) (2 * dxU1 / domainWidth) * (0.5 - x/domainWidth) .* z;
    
    dxXiCont = @(t,x) 2 * dxXi/domainWidth * (x/domainWidth - 0.5);
    dtHCont = @(t,x) dtXi * ones(size(x));
    
    dtU1Cont = @(t,x,z) zeros(size(x));
    dxU1Cont = @(t,x,z) 2 * dxU1 / domainWidth * (x/domainWidth - 0.5);
    dzU1Cont = @(t,x,z) zeros(size(x));
    dzU2Cont = @(t,x,z) 2 * dxU1 / domainWidth * (0.5 - x/domainWidth);
    
    dxdxU1Cont = @(t,x,z) zeros(size(x));
    dxdzU1Cont = @(t,x,z) zeros(size(x));
    dzdzU1Cont = @(t,x,z) zeros(size(x));
    
    dxU1zIntCont = @(t,x) dxU1Cont(t,x,0) .* hCont(t,x) + u1Cont(t,x,0) .* (dxXiCont(t,x) - dxZb);
    
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
      @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    dxzDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
      @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
    
  case 'quadratic' %
    domainWidth = 100;
    idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idBdrRiem = [2,4];
%     idLand = -1; idOS = -1; idRiv = -1; idRad = -1; idRiem = -1;
    %     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = -1;
    
    gConst = 10;
    xi0Cont = @(x) zeros(size(x));
    dtXi = 0.01;
    dxZb = -0.005;
    dxXi = -0.3;
    dxU1 = 0.01;
    dzU2 = -dxU1;
    rho = 0.1;
    
    xiCont = @(t,x) dxXi * (x/domainWidth - 0.5).^2 + dtXi * t;
    zBotCont = @(x) -2 + dxZb * x;
    
    hCont = @(t,x) xiCont(t,x) - zBotCont(x);
    u1Cont = @(t,x,z) dxU1 * x;
    u2Cont = @(t,x,z) dzU2 * z;
    
    dxXiCont = @(t,x) 2 * dxXi/domainWidth * (x/domainWidth - 0.5);
    dtHCont = @(t,x) dtXi * ones(size(x));
    
    dtU1Cont = @(t,x,z) zeros(size(x));
    dxU1Cont = @(t,x,z) dxU1 * ones(size(x));
    dzU1Cont = @(t,x,z) zeros(size(x));
    dzU2Cont = @(t,x,z) dzU2 * ones(size(x));
    
    dxdxU1Cont = @(t,x,z) zeros(size(x));
    dxdzU1Cont = @(t,x,z) zeros(size(x));
    dzdzU1Cont = @(t,x,z) zeros(size(x));
    
    u1zIntCont = @(t,x) dxU1 * x .* hCont(t,x);
    dxU1zIntCont = @(t,x) dxU1 * hCont(t,x) + dxU1 * x .* (dxXiCont(t,x) - dxZb);
    
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
      @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    dxzDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
      @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
    
  case 'test_h'
    domainWidth = 100;
    idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idBdrRiem = [2,4];
    %     idLand = -1; idOS = -1; idRiv = -1; idRad = -1; idRiem = -1;
    %     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = -1;
    
    gConst = 10;
    xi0Cont = @(x) zeros(size(x));
    dxZb = -0.005;
    omega = 0.1;
    eta = 0.1;
    epsilon = 0.01;
    rho = 0;
    
    xiCont = @(t,x) epsilon * sin(omega * x + eta * t);
    zBotCont = @(x) -2 + dxZb * x;
    
    hCont = @(t,x) xiCont(t,x) - zBotCont(x);
    u1Cont = @(t,x,z) 0.1 * ones(size(x));
    u2Cont = @(t,x,z) zeros(size(x));
    
    dxXiCont = @(t,x) epsilon * omega * cos(omega * x + eta * t);
    dtHCont = @(t,x) epsilon * eta * cos(omega * x + eta * t);
    
    dtU1Cont = @(t,x,z) zeros(size(x));
    dxU1Cont = @(t,x,z) zeros(size(x));
    dzU1Cont = @(t,x,z) zeros(size(x));
    dzU2Cont = @(t,x,z) zeros(size(x));
    
    dxdxU1Cont = @(t,x,z) zeros(size(x));
    dxdzU1Cont = @(t,x,z) zeros(size(x));
    dzdzU1Cont = @(t,x,z) zeros(size(x));
    
    dxU1zIntCont = @(t,x) 0.1 * (dxXiCont(t,x) - dxZb);
    
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
      @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    dxzDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
      @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
    
  case 'test_u'
    domainWidth = 100;
    %     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = [2,4];
    %     idLand = -1; idOS = -1; idRiv = -1; idRad = -1; idRiem = -1;
    %     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = -1;
    idLand = -1; idOS = -1; idRiv = 4; idRad = -1; idBdrRiem = -1;
    
    gConst = 10;
    xi0Cont = @(x) zeros(size(x));
    dtXi = 0.025;
    dxZb = 0;%-0.005;
    omega = 0.1;
    delta = 0.1;
    rho = 0;
    
    xiCont = @(t,x) dtXi * t * ones(size(x));
    zBotCont = @(x) -2 + dxZb * x;
    
    hCont = @(t,x) xiCont(t,x) - zBotCont(x);
    u1Cont = @(t,x,z) delta * sin(omega * t + 4 * z);
    u2Cont = @(t,x,z) 0.1 * ones(size(x));
    
    dxXiCont = @(t,x) zeros(size(x));
    dtHCont = @(t,x) dtXi * ones(size(x));
    
    dtU1Cont = @(t,x,z) delta * omega * cos(omega * t + 4 * z);
    dxU1Cont = @(t,x,z) zeros(size(x));
    dzU1Cont = @(t,x,z) delta * 4 * cos(omega * t + 4 * z);
    dzU2Cont = @(t,x,z) zeros(size(x));
    
    dxdxU1Cont = @(t,x,z) zeros(size(x));
    dxdzU1Cont = @(t,x,z) zeros(size(x));
    dzdzU1Cont = @(t,x,z) zeros(size(x));
    
    u1zIntCont = @(t,x) 0.25 * delta * (-cos(omega * t + 4 * xiCont(t,x)) + cos(omega * t + 4 * zBotCont(x)));
    dxU1zIntCont = @(t,x) -delta * sin(omega * t + zBotCont(x)) * dxZb;
    
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
      @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    dxzDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
      @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
    
  case 'test_w'
    domainWidth = 100;
    %     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = [2,4];
    %     idLand = -1; idOS = -1; idRiv = -1; idRad = -1; idRiem = -1;
    %     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = -1;
    idLand = -1; idOS = -1; idRiv = 4; idRad = -1; idBdrRiem = -1;
    
    gConst = 10;
    xi0Cont = @(x) zeros(size(x));
    dtXi = 0.025;
    dxZb = -0.005;
    omega = 0.1;
    delta = 0.1;
    rho = 0.001;
    
    xiCont = @(t,x) dtXi * t * ones(size(x));
    zBotCont = @(x) -2 + dxZb * x;
    
    hCont = @(t,x) xiCont(t,x) - zBotCont(x);
    u1Cont = @(t,x,z) 0.1 * ones(size(x));
    u2Cont = @(t,x,z) delta * sin(omega * t + 4 * x);
    
    dxXiCont = @(t,x) zeros(size(x));
    dtHCont = @(t,x) dtXi * ones(size(x));
    
    dtU1Cont = @(t,x,z) zeros(size(x));
    dxU1Cont = @(t,x,z) zeros(size(x));
    dzU1Cont = @(t,x,z) zeros(size(x));
    dzU2Cont = @(t,x,z) zeros(size(x));
    
    dxdxU1Cont = @(t,x,z) zeros(size(x));
    dxdzU1Cont = @(t,x,z) zeros(size(x));
    dzdzU1Cont = @(t,x,z) zeros(size(x));
    
    u1zIntCont = @(t,x) 0.1 * hCont(t,x);
    dxU1zIntCont = @(t,x) zeros(size(x));
    
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
      @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    dxzDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
      @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
    
    
  case 'test_uw'
    domainWidth = 100;
    %     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = [2,4];
    idLand = -1; idOS = -1; idRiv = -1; idRad = -1; idBdrRiem = -1;
    %     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = [2,4]; idRiem = -1;
    %     idLand = -1; idOS = -1; idRiv = [2,4]; idRad = -1; idRiem = -1;
    
    gConst = 10;
    xi0Cont = @(x) zeros(size(x));
    dtXi = 0.0;
    dxZb = -0.005;
    omega = 0.1;
    theta = 0.1;
    delta = 0.1;
    rho = 0.0;%0.001;
    
    xiCont = @(t,x) dtXi * t * ones(size(x));
    zBotCont = @(x) -2 + dxZb * x;
    
    hCont = @(t,x) xiCont(t,x) - zBotCont(x);
    u1Cont = @(t,x,z) delta * (z - zBotCont(x)) .* sin(omega * x + theta * t);
    u2Cont = @(t,x,z) delta * dxZb * z .* sin(omega * x + theta * t) - 0.5 * delta * omega * (z - zBotCont(x)).^2 .* cos(omega * x + theta * t);
    
    dxXiCont = @(t,x) zeros(size(x));
    dtHCont = @(t,x) dtXi * ones(size(x));
    
    dtU1Cont = @(t,x,z) delta * theta * (z - zBotCont(x)) .* cos(omega * x + theta * t);
    dxU1Cont = @(t,x,z) -delta * dxZb * sin(omega * x + theta * t) + delta * omega * (z - zBotCont(x)) .* cos(omega * x + theta * t);
    dzU1Cont = @(t,x,z) delta * sin(omega * x + theta * t);
    dzU2Cont = @(t,x,z) delta * dxZb * sin(omega * x + theta * t) - delta * omega * (z - zBotCont(x)) .* cos(omega * x + theta * t);
    
    dxdxU1Cont = @(t,x,z) -delta * omega * ( 2 * dxZb * cos(omega * x + theta * t) + omega * (z - zBotCont(x)) .* sin(omega * x + theta * t) );
    dxdzU1Cont = @(t,x,z) delta * omega * cos(omega * x + theta * t);
    dzdzU1Cont = @(t,x,z) zeros(size(x));
    
    u1zIntCont = @(t,x) 0.5 * delta * sin(omega * x + theta * t) .* (xiCont(t,x) - zBotCont(x)).^2;
    dxU1zIntCont = @(t,x) 0.5 * delta * omega * cos(omega * x + theta * t) .* (xiCont(t,x) - zBotCont(x)).^2 + ...
      delta * sin(omega * x + theta * t) .* (xiCont(t,x) - zBotCont(x)) .* (dxXiCont(t,x) - dxZb);
    
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
      @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    dxzDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
      @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
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
  fuCont = @(t,x,z) dtU1Cont(t,x,z) + 2 * u1Cont(t,x,z) .* dxU1Cont(t,x,z) + ...
    dzU1Cont(t,x,z) .* u2Cont(t,x,z) + u1Cont(t,x,z) .* dzU2Cont(t,x,z) - ( ...
    DCont{1,1}(t,x,z) .* dxdxU1Cont(t,x,z) + dxzDCont{1,1}(t,x,z) .* dxU1Cont(t,x,z) + ...
    DCont{1,2}(t,x,z) .* dxdzU1Cont(t,x,z) + dxzDCont{1,2}(t,x,z) .* dzU1Cont(t,x,z) + ...
    DCont{2,1}(t,x,z) .* dxdzU1Cont(t,x,z) + dxzDCont{2,1}(t,x,z) .* dxU1Cont(t,x,z) + ...
    DCont{2,2}(t,x,z) .* dzdzU1Cont(t,x,z) + dxzDCont{2,2}(t,x,z) .* dzU1Cont(t,x,z) ) + ...
    gConst * dxXiCont(t,x);

  problemData = setdefault(problemData, 'fhCont', fhCont);
  problemData = setdefault(problemData, 'fuCont', fuCont);

  problemData = setdefault(problemData, 'h0Cont', @(x1) hCont(t0, x1));
  problemData = setdefault(problemData, 'u10Cont', @(x1,x2) u1Cont(t0, x1, x2));

  problemData = setdefault(problemData, 'hDCont', hCont);
  problemData = setdefault(problemData, 'u1DCont', u1Cont);
  problemData = setdefault(problemData, 'u2DCont', u2Cont);
  problemData = setdefault(problemData, 'q1DCont', @(t,x,z) -DCont{1,1}(t,x,z) .* dxU1Cont(t,x,z) - DCont{1,2}(t,x,z) .* dzU1Cont(t,x,z));
  problemData = setdefault(problemData, 'q2DCont', @(t,x,z) -DCont{2,1}(t,x,z) .* dxU1Cont(t,x,z) - DCont{2,2}(t,x,z) .* dzU1Cont(t,x,z));
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
