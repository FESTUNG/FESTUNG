function problemData = getTestcase(problemData, problemName)
switch problemName  
  case 'convergence'
    domainWidth = 100;
    idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; 
%     idLand = -1; idOS = -1; idRiv = -1; idRad = -1; 
    
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
    
    u1hCont = @(t,x) 0.5 * delta * sin(omega * (x+t)) .* (xiCont(t,x) - zBotCont(x)).^2;
    dxU1hCont = @(t,x) 0.5 * delta * omega * cos(omega * (x+t)) .* (xiCont(t,x) - zBotCont(x)).^2 + ...
                        delta * sin(omega * (x+t)) .* (xiCont(t,x) - zBotCont(x)) .* (dxXiCont(t,x) - dxZb);
                        
    DCont = { @(t,x,z) rho * ones(size(x)), @(t,x,z) zeros(size(x)); ...
              @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
    dxzDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
                 @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
  
  case 'coupling'
    domainWidth = 100;
    idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1;
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
    idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1;
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
        
  case 'test'
    domainWidth = 100;
    idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; 
%     idLand = -1; idOS = -1; idRiv = -1; idRad = -1; 
    
    gConst = 10;
    xi0Cont = @(x) zeros(size(x));
    dxZb = 0;
    
    xiCont = @(t,x) 0.1 * sin(0.01*(x+t));
    zBotCont = @(x) -2 + dxZb * x;
    
    hCont = @(t,x) xiCont(t,x) - zBotCont(x);
    u1Cont = @(t,x,z) ones(size(x)); %sin(z+0.01*t);
    u2Cont = @(t,x,z) zeros(size(x));
    
    dxXiCont = @(t,x) 0.1*0.01*cos(0.01*(x+t));
    dtHCont = @(t,x) 0.1*0.01*cos(0.01*(x+t));
    
    dtU1Cont = @(t,x,z) zeros(size(x));%0.01*cos(z+0.01*t);
    dxU1Cont = @(t,x,z) zeros(size(x));
    dzU1Cont = @(t,x,z) zeros(size(x));%cos(z+0.01*t);
    dzU2Cont = @(t,x,z) zeros(size(x));
    
    dxdxU1Cont = @(t,x,z) zeros(size(x));
    dxdzU1Cont = @(t,x,z) zeros(size(x));
    dzdzU1Cont = @(t,x,z) zeros(size(x));%sin(z+0.01*t);
    
    u1hCont = @(t,x) hCont(t,x);%cos(xiCont(t,x)+0.01*t) - cos(zBotCont(x)+0.01*t);
    dxU1hCont = @(t,x) dxXiCont(t,x);%-sin(xiCont(t,x)+0.01*t).*dxXiCont(t,x);

    DCont = cellfun(@(c) @(t,x,z) c * ones(size(x)), {0.01, 0; 0, 0.01}, 'UniformOutput', false);
    dxzDCont = cellfun(@(c) @(t,x,z) c * ones(size(x)), {0, 0; 0, 0}, 'UniformOutput', false);

end % switch

problemData = setdefault(problemData, 'domainWidth', domainWidth);
problemData = setdefault(problemData, 'xi0Cont', xi0Cont);
problemData = setdefault(problemData, 'zBotCont', zBotCont);
problemData = setdefault(problemData, 'idLand', idLand);
problemData = setdefault(problemData, 'idOS', idOS);
problemData = setdefault(problemData, 'idRiv', idRiv);
problemData = setdefault(problemData, 'idRad', idRad);

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
