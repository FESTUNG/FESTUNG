function problemData = getTestcase(problemData, problemName)
switch problemName  
  case 'convergence'
    domainWidth = 100;
%     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = [2,4];
%     idLand = -1; idOS = -1; idRiv = -1; idRad = -1; idRiem = -1;
    idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = [2,4]; idRiem = -1;
    
    gConst = 10;
    xi0Cont = @(x) zeros(size(x));
    dxZb = -0.005;
    omega = 0.1;
    delta = 0.1;
    epsilon = 0.01;
    rho = 0.001;
    
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
    
    u1hCont = @(t,x) 0.5 * delta * sin(omega * (x+t)) .* (xiCont(t,x).^2 - zBotCont(x).^2) - delta * zBotCont(x) .* hCont(t,x) .* sin(omega * (x+t));
    dxU1hCont = @(t,x) 0.5 * delta * omega * cos(omega * (x+t)) .* (xiCont(t,x) - zBotCont(x)).^2 + ...
                        delta * sin(omega * (x+t)) .* (xiCont(t,x) - zBotCont(x)) .* (dxXiCont(t,x) - dxZb);
                        
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
              @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    dxzDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
                 @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
  
  case 'coupling'
    domainWidth = 100;
    idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = [2,4];
%     idLand = -1; idOS = -1; idRiv = -1; idRad = -1;
    
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
    
    u1hCont = @(t,x) -k * a * d / b * sin(b*x+c*t) .* ( sin(d*xiCont(t,x)) - sin(d*zBotCont(x)) ) + 0.5 * f * ( xiCont(t,x).^2 - zBotCont(x).^2 );
    dxU1hCont = @(t,x) -k * a * d * cos(b*x+c*t) .* ( sin(d*xiCont(t,x)) - sin(d*zBotCont(x)) ) - ...
                        k * a * d^2 / b * sin(b*x+c*t) .* ( cos(d*xiCont(t,x)) .* dxXiCont(t,x) - cos(d*zBotCont(x)) * dxZb ) + ...
                        f * ( xiCont(t,x) .* dxXiCont(t,x) - zBotCont(x) * dxZb );
                        
    DCont = cellfun(@(c) @(t,x,z) c * ones(size(x)), {0.01, 0; 0, 0.01}, 'UniformOutput', false);                        
    dxzDCont = cellfun(@(c) @(t,x,z) c * ones(size(x)), {0, 0; 0, 0}, 'UniformOutput', false);
    
  case 'utbest'
    domainWidth = 100;
    idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = [2,4];
%     idLand = -1; idOS = -1; idRiv = -1; idRad = -1; 
    
    gConst = 10;
    xi0Cont = @(x) zeros(size(x));
    dxZb = 0.005;
    omega = 0.01;
    d = 0.1;
    e = 0.01;
    D = 0.1;
    
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
    
    u1hCont = @(t,x) 0.5 * d * sin(omega * (x+t)) .* (xiCont(t,x) - zBotCont(x)).^2;
    dxU1hCont = @(t,x) 0.5 * d * omega * cos(omega * (x+t)) .* (xiCont(t,x) - zBotCont(x)).^2 + ...
                        d * sin(omega * (x+t)) .* (xiCont(t,x) - zBotCont(x)) .* (dxXiCont(t,x) - dxZb);
                        
    DCont = cellfun(@(c) @(t,x,z) c * ones(size(x)), {D, 0; 0, D}, 'UniformOutput', false);
    dxzDCont = cellfun(@(c) @(t,x,z) c * ones(size(x)), {0, 0; 0, 0}, 'UniformOutput', false);
        
  case 'utbest_sinus'
    domainWidth = 100;
    idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = [2,4];
%     idLand = -1; idOS = -1; idRiv = -1; idRad = -1; idRiem = -1;
    
    gConst = 10;
    xi0Cont = @(x) zeros(size(x));
    dxZb = 0;
    omega = 0.01;
    d = 0.1;
    e = 0.01;
    t_coef = 1;
    D = 0.1;
    
    xiCont = @(t,x) e * sin(omega * (x+t_coef*t));
    zBotCont = @(x) -2 + dxZb * x;
    
    hCont = @(t,x) xiCont(t,x) - zBotCont(x);
    u1Cont = @(t,x,z) d * (sin(omega * (x+t_coef*t)) + sin(omega * (z+t_coef*t)));
    u2Cont = @(t,x,z) -d * omega * z .* cos(omega * (x+t_coef*t));
    
    dxXiCont = @(t,x) e * omega * cos(omega * (x+t_coef*t));
    dtHCont = @(t,x)  e * omega * t_coef * cos(omega * (x+t_coef*t));
    
    dtU1Cont = @(t,x,z) d * omega * t_coef * (cos(omega * (x+t_coef*t)) + cos(omega * (z+t_coef*t)));
    dxU1Cont = @(t,x,z) d * omega * cos(omega * (x+t_coef*t));
    dzU1Cont = @(t,x,z) d * omega * cos(omega * (z+t_coef*t));
    dzU2Cont = @(t,x,z) -d * omega * cos(omega * (x+t_coef*t));
    
    dxdxU1Cont = @(t,x,z) -d * omega^2 * sin(omega * (x+t_coef*t));
    dxdzU1Cont = @(t,x,z) zeros(size(x));
    dzdzU1Cont = @(t,x,z) zeros(size(x));
    
    u1hCont = @(t,x) d * sin(omega * (x+t_coef*t)) .* hCont(t,x) - d/omega * (cos(omega * (xiCont(t,x)+t_coef*t)) - cos(omega * (zBotCont(x)+t_coef*t)));
    dxU1hCont = @(t,x) d * omega * cos(omega * (x+t_coef*t)) .* hCont(t,x) + d * sin(omega * (x+t_coef*t)) .* dxXiCont(t,x) + ...
                      d * sin(omega * (xiCont(t,x)+t_coef*t)) .* dxXiCont(t,x);
                        
    DCont = cellfun(@(c) @(t,x,z) c * ones(size(x)), {D, 0; 0, D}, 'UniformOutput', false);
    dxzDCont = cellfun(@(c) @(t,x,z) c * ones(size(x)), {0, 0; 0, 0}, 'UniformOutput', false);
        
  case 'test_no_diffusion'
    domainWidth = 100;
%     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = [2,4]; idRiem = [2,4];
    idLand = -1; idOS = -1; idRiv = -1; idRad = -1; idRiem = -1;
    
    gConst = 10;
    xi0Cont = @(x) zeros(size(x));
    dxZb = 0.005;%0;
    
    xiCont = @(t,x) 0.01 * sin(0.01*(100-x));%0.005 * x; %0.01 * sin(0.1*(x));
    zBotCont = @(x) -2 + dxZb * x;
    
    hCont = @(t,x) xiCont(t,x) - zBotCont(x);
    u1Cont = @(t,x,z) sin(z+0.01*t);%zeros(size(x)); %sin(z+0.01*t);
    u2Cont = @(t,x,z) zeros(size(x));
    
    dxXiCont = @(t,x) -0.01*0.01*cos(0.01*(100-x));%0.005 * ones(size(x)); %0.1*0.01*cos(0.1*(x));
    dtHCont = @(t,x) zeros(size(x));
    
    dtU1Cont = @(t,x,z) 0.01*cos(z+0.01*t);%zeros(size(x));%0.01*cos(z+0.01*t);
    dxU1Cont = @(t,x,z) zeros(size(x));
    dzU1Cont = @(t,x,z) cos(z+0.01*t);%zeros(size(x));%cos(z+0.01*t);
    dzU2Cont = @(t,x,z) zeros(size(x));
    
    dxdxU1Cont = @(t,x,z) zeros(size(x));
    dxdzU1Cont = @(t,x,z) zeros(size(x));
    dzdzU1Cont = @(t,x,z) -sin(z+0.01*t);%zeros(size(x));%-sin(z+0.01*t);
    
    u1hCont = @(t,x) -cos(xiCont(t,x)+0.01*t) + cos(zBotCont(x)+0.01*t);%zeros(size(x));%cos(xiCont(t,x)+0.01*t) - cos(zBotCont(x)+0.01*t);
    dxU1hCont = @(t,x) sin(xiCont(t,x)+0.01*t).*dxXiCont(t,x);%zeros(size(x));%-sin(xiCont(t,x)+0.01*t).*dxXiCont(t,x);

    DCont = cellfun(@(c) @(t,x,z) c * ones(size(x)), {0.0, 0; 0, 0.0}, 'UniformOutput', false);
    dxzDCont = cellfun(@(c) @(t,x,z) c * ones(size(x)), {0, 0; 0, 0}, 'UniformOutput', false);

  case 'constant' % OK (ausser mit Riemann-Loeser am Rand)
    domainWidth = 100;
%     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = [2,4];
%     idLand = -1; idOS = -1; idRiv = -1; idRad = -1; idRiem = -1;
    idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = -1;
    
    gConst = 10;
    xi0Cont = @(x) zeros(size(x));
    dtXi = 0.01;
    rho = 0;
    
    xiCont = @(t,x) dtXi * t * ones(size(x));
    zBotCont = @(x) -2 * ones(size(x));
    
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
    
    u1hCont = @(t,x) 0.1 * hCont(t,x);
    dxU1hCont = @(t,x) zeros(size(x));
                        
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
              @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    dxzDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
                 @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
    
  case 'linear_h' % OK (ausser mit Riemann-Loeser am Rand)
    domainWidth = 100;
    idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = [2,4];
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
    
    u1hCont = @(t,x) 0.1 * hCont(t,x);
    dxU1hCont = @(t,x) 0.1 * (dxXiCont(t,x) - dxZb);
                        
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
              @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    dxzDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
                 @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
    
  case 'z-linear_u' % OK (ausser mit Riemann-Loeser am Rand)
    domainWidth = 100;
%     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = [2,4];
%     idLand = -1; idOS = -1; idRiv = -1; idRad = -1; idRiem = -1;
    idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = -1;
    
    gConst = 10;
    xi0Cont = @(x) zeros(size(x));
    dzU1 = 0.1;
    rho = 0;
    
    xiCont = @(t,x) zeros(size(x));
    zBotCont = @(x) -2 * ones(size(x));
    
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
    
    u1hCont = @(t,x) -0.2 * ones(size(x));
    dxU1hCont = @(t,x) zeros(size(x));
                        
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
              @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    dxzDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
                 @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
    
  case 'x-linear_w' % OK (ausser mit Riemann-Loeser am Rand)
    domainWidth = 100;
%     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = [2,4];
%     idLand = -1; idOS = -1; idRiv = -1; idRad = -1; idRiem = -1;
    idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = -1;
    
    gConst = 10;
    xi0Cont = @(x) zeros(size(x));
    dxU2 = 0.01;
    rho = 0;
    
    xiCont = @(t,x) zeros(size(x));
    zBotCont = @(x) -2 * ones(size(x));
    
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
    
    u1hCont = @(t,x) 0.2 * ones(size(x));
    dxU1hCont = @(t,x) zeros(size(x));
                        
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
              @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    dxzDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
                 @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
    
  case 'linear_vel' % OK (ausser mit Riemann-Loeser am Rand)
    domainWidth = 100;
%     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = [2,4];
%     idLand = -1; idOS = -1; idRiv = -1; idRad = -1; idRiem = -1;
    idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = -1;
    
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
    
    u1hCont = @(t,x) 2 * dxU1 * x;
    dxU1hCont = @(t,x) 2 * dxU1 * ones(size(x));
                        
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
              @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    dxzDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
                 @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
    
  case 'linear' % OK (ausser mit Riemann-Loeser am Rand)
    domainWidth = 100;
%     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = [2,4];
%     idLand = -1; idOS = -1; idRiv = -1; idRad = -1; idRiem = -1;
%     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = -1;
    idLand = -1; idOS = -1; idRiv = [2,4]; idRad = -1; idRiem = [2,4];
    
    gConst = 10;
    xi0Cont = @(x) zeros(size(x));
    dtXi = 0;%0.01;
    dxZb = 0;%-0.005;
    dxXi = 0;%0.005;
    dxU1 = 0.01;
    dzU2 = -dxU1;
    rho = 0;
    
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
    
    u1hCont = @(t,x) dxU1 * x .* hCont(t,x);
    dxU1hCont = @(t,x) dxU1 * hCont(t,x) + dxU1 * x .* (dxXiCont(t,x) - dxZb);
                        
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
              @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    dxzDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
                 @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
    
  case 'quadratic_h' % 
    domainWidth = 100;
%     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = [2,4];
    idLand = -1; idOS = -1; idRiv = -1; idRad = -1; idRiem = -1;
%     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = -1;
%     idLand = -1; idOS = -1; idRiv = [2,4]; idRad = -1; idRiem = -1;
    
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
    
    u1hCont = @(t,x) dxU1 * x .* hCont(t,x);
    dxU1hCont = @(t,x) dxU1 * hCont(t,x) + dxU1 * x .* (dxXiCont(t,x) - dxZb);
                        
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
              @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    dxzDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
                 @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
    
  case 'quadratic_u' % 
    domainWidth = 100;
%     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = [2,4];
%     idLand = -1; idOS = -1; idRiv = -1; idRad = -1; idRiem = -1;
%     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = -1;
    idLand = -1; idOS = -1; idRiv = 4; idRad = 2; idRiem = -1;
    
    gConst = 10;
    xi0Cont = @(x) zeros(size(x));
    dtXi = 0;%0.01;
    dxZb = 0;%-0.005;
    dxXi = 0;%-0.3;
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
    
    u1hCont = @(t,x) dxU1 * x .* hCont(t,x);
    dxU1hCont = @(t,x) dxU1 * hCont(t,x) + dxU1 * x .* (dxXiCont(t,x) - dxZb);
                        
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
              @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    dxzDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
                 @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
    
  case 'quadratic' % 
    domainWidth = 100;
%     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = [2,4];
    idLand = -1; idOS = -1; idRiv = -1; idRad = -1; idRiem = -1;
%     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = -1;
%     idLand = -1; idOS = -1; idRiv = [2,4]; idRad = -1; idRiem = -1;
    
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
    
    u1hCont = @(t,x) dxU1 * x .* hCont(t,x);
    dxU1hCont = @(t,x) dxU1 * hCont(t,x) + dxU1 * x .* (dxXiCont(t,x) - dxZb);
                        
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
              @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    dxzDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
                 @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
    
  case 'test_h'
    domainWidth = 100;
%     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = [2,4];
    idLand = -1; idOS = -1; idRiv = -1; idRad = -1; idRiem = -1;
%     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = -1;
    
    gConst = 10;
    xi0Cont = @(x) zeros(size(x));
    dxZb = -0.005;
    omega = 0.1;
    epsilon = 0.01;
    rho = 0.001;
    
    xiCont = @(t,x) epsilon * sin(omega * (x+t));
    zBotCont = @(x) -2 + dxZb * x;
    
    hCont = @(t,x) xiCont(t,x) - zBotCont(x);
    u1Cont = @(t,x,z) 0.1 * ones(size(x));
    u2Cont = @(t,x,z) zeros(size(x));
    
    dxXiCont = @(t,x) epsilon * omega * cos(omega * (x+t));
    dtHCont = @(t,x) epsilon * omega * cos(omega * (x+t));
        
    dtU1Cont = @(t,x,z) zeros(size(x));
    dxU1Cont = @(t,x,z) zeros(size(x));
    dzU1Cont = @(t,x,z) zeros(size(x));
    dzU2Cont = @(t,x,z) zeros(size(x));
    
    dxdxU1Cont = @(t,x,z) zeros(size(x));
    dxdzU1Cont = @(t,x,z) zeros(size(x));
    dzdzU1Cont = @(t,x,z) zeros(size(x));
    
    u1hCont = @(t,x) 0.1 * hCont(t,x);
    dxU1hCont = @(t,x) 0.1 * (dxXiCont(t,x) - dxZb);
                        
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
              @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    dxzDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
                 @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
    
  case 'test_u'
    domainWidth = 100;
%     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = [2,4];
%     idLand = -1; idOS = -1; idRiv = -1; idRad = -1; idRiem = -1;
%     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = -1;
    idLand = -1; idOS = -1; idRiv = 4; idRad = -1; idRiem = -1;
    
    gConst = 10;
    xi0Cont = @(x) zeros(size(x));
    dtXi = 0.025;
    dxZb = 0;%-0.005;
    omega = 0.1;
    delta = 0.1;
    rho = 0.001;
    
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
    
    u1hCont = @(t,x) 0.25 * delta * (-cos(omega * t + 4 * xiCont(t,x)) + cos(omega * t + 4 * zBotCont(x)));
    dxU1hCont = @(t,x) -delta * sin(omega * t + zBotCont(x)) * dxZb;
                        
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
              @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    dxzDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
                 @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
               
  case 'test_w'
    domainWidth = 100;
%     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = [2,4];
%     idLand = -1; idOS = -1; idRiv = -1; idRad = -1; idRiem = -1;
%     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = -1;
    idLand = -1; idOS = -1; idRiv = 4; idRad = -1; idRiem = -1;
    
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
    
    u1hCont = @(t,x) 0.1 * hCont(t,x);
    dxU1hCont = @(t,x) zeros(size(x));
                        
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
              @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    dxzDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
                 @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
    

  case 'test_uw'
    domainWidth = 100;
%     idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idRiem = [2,4];
    idLand = -1; idOS = -1; idRiv = -1; idRad = -1; idRiem = -1;
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
    
    u1hCont = @(t,x) 0.5 * delta * sin(omega * x + theta * t) .* (xiCont(t,x) - zBotCont(x)).^2;
    dxU1hCont = @(t,x) 0.5 * delta * omega * cos(omega * x + theta * t) .* (xiCont(t,x) - zBotCont(x)).^2 + ...
                        delta * sin(omega * x + theta * t) .* (xiCont(t,x) - zBotCont(x)) .* (dxXiCont(t,x) - dxZb);
    
    DCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
              @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    dxzDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
                 @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
end % switch

problemData = setdefault(problemData, 'domainWidth', domainWidth);
problemData = setdefault(problemData, 'xi0Cont', xi0Cont);
problemData = setdefault(problemData, 'zBotCont', zBotCont);
problemData = setdefault(problemData, 'idLand', idLand);
problemData = setdefault(problemData, 'idOS', idOS);
problemData = setdefault(problemData, 'idRiv', idRiv);
problemData = setdefault(problemData, 'idRad', idRad);
problemData = setdefault(problemData, 'idRiem', idRiem);

problemData = setdefault(problemData, 'gConst', gConst);
problemData = setdefault(problemData, 'hCont', hCont);
problemData = setdefault(problemData, 'u1Cont', u1Cont);
problemData = setdefault(problemData, 'u2Cont', u2Cont);
problemData = setdefault(problemData, 'DCont', DCont);

fhCont = @(t,x) dtHCont(t,x) + dxU1hCont(t,x);
fuCont = @(t,x,z) dtU1Cont(t,x,z) + 2 * u1Cont(t,x,z) .* dxU1Cont(t,x,z) + ...
                  dzU1Cont(t,x,z) .* u2Cont(t,x,z) + u1Cont(t,x,z) .* dzU2Cont(t,x,z) - ( ...
                    DCont{1,1}(t,x,z) .* dxdxU1Cont(t,x,z) + dxzDCont{1,1}(t,x,z) .* dxU1Cont(t,x,z) + ...
                    DCont{1,2}(t,x,z) .* dxdzU1Cont(t,x,z) + dxzDCont{1,2}(t,x,z) .* dzU1Cont(t,x,z) + ...
                    DCont{2,1}(t,x,z) .* dxdzU1Cont(t,x,z) + dxzDCont{2,1}(t,x,z) .* dxU1Cont(t,x,z) + ...
                    DCont{2,2}(t,x,z) .* dzdzU1Cont(t,x,z) + dxzDCont{2,2}(t,x,z) .* dzU1Cont(t,x,z) ) + ...
                    gConst * dxXiCont(t,x);

problemData = setdefault(problemData, 'fhCont', fhCont);
problemData = setdefault(problemData, 'fuCont', fuCont);

problemData = setdefault(problemData, 'hDCont', hCont);
problemData = setdefault(problemData, 'u1DCont', u1Cont);
problemData = setdefault(problemData, 'u2DCont', u2Cont);
problemData = setdefault(problemData, 'q1DCont', @(t,x,z) -DCont{1,1}(t,x,z) .* dxU1Cont(t,x,z) - DCont{1,2}(t,x,z) .* dzU1Cont(t,x,z));
problemData = setdefault(problemData, 'q2DCont', @(t,x,z) -DCont{2,1}(t,x,z) .* dxU1Cont(t,x,z) - DCont{2,2}(t,x,z) .* dzU1Cont(t,x,z));
problemData = setdefault(problemData, 'uhDCont', @(t,x) u1hCont(t,x));
end % function
