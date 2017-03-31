function [problemData, domainWidth, xi0Cont, zBotCont, idLand, idOS, idRiv, idRad, idCoupling] = getTestcase(problemData, problemName)
switch problemName  
  case 'convergence'
    domainWidth = 100;
    idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idCoupling = -1;
%     idLand = -1; idOS = -1; idRiv = -1; idRad = -1; idCoupling = -1;
    
    problemData.gConst = 10;
    xi0Cont = @(x) zeros(size(x));
    dxZb = -0.005;
    omega = 0.1;
    delta = 0.1;
    epsilon = 0.02;
    rho = 0.001;
    
    xiCont = @(t,x) epsilon * sin(omega * (x+t));
    zBotCont = @(x) -2 + dxZb * x;
    
    problemData.hCont = @(t,x) xiCont(t,x) - zBotCont(x);
    problemData.u1Cont = @(t,x,z) delta * (z - zBotCont(x)) .* sin(omega * (x+t));
    problemData.u2Cont = @(t,x,z) delta * dxZb * z .* sin(omega * (x+t)) - 0.5 * delta * omega * (z - zBotCont(x)).^2 .* cos(omega * (x+t));
    
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
                        
    problemData.DCont = { @(t,x,z) rho * ones(size(x)), @(t,x,z) zeros(size(x)); ...
                          @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
                        
    dxzDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
                 @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
    
    problemData.fhCont = @(t,x) dtHCont(t,x) + dxU1hCont(t,x);
    problemData.fuCont = @(t,x,z) dtU1Cont(t,x,z) + 2 * problemData.u1Cont(t,x,z) .* dxU1Cont(t,x,z) + ...
                            dzU1Cont(t,x,z) .* problemData.u2Cont(t,x,z) + problemData.u1Cont(t,x,z) .* dzU2Cont(t,x,z) - ( ...
                              problemData.DCont{1,1}(t,x,z) .* dxdxU1Cont(t,x,z) + dxzDCont{1,1}(t,x,z) .* dxU1Cont(t,x,z) + ...
                              problemData.DCont{1,2}(t,x,z) .* dxdzU1Cont(t,x,z) + dxzDCont{1,2}(t,x,z) .* dzU1Cont(t,x,z) + ...
                              problemData.DCont{2,1}(t,x,z) .* dxdzU1Cont(t,x,z) + dxzDCont{2,1}(t,x,z) .* dxU1Cont(t,x,z) + ...
                              problemData.DCont{2,2}(t,x,z) .* dzdzU1Cont(t,x,z) + dxzDCont{2,2}(t,x,z) .* dzU1Cont(t,x,z) ) + ...
                              problemData.gConst * dxXiCont(t,x);
    
    problemData.hDCont = problemData.hCont;
    problemData.u1DCont = problemData.u1Cont;
    problemData.u2DCont = problemData.u2Cont;
    problemData.q1DCont = @(t,x,z) -problemData.DCont{1,1}(t,x,z) .* dxU1Cont(t,x,z) - problemData.DCont{1,2}(t,x,z) .* dzU1Cont(t,x,z);
    problemData.q2DCont = @(t,x,z) -problemData.DCont{2,1}(t,x,z) .* dxU1Cont(t,x,z) - problemData.DCont{2,2}(t,x,z) .* dzU1Cont(t,x,z);
    problemData.uhDCont = @(t,x) u1hCont(t,x);
  
  case 'convergence2'
    domainWidth = 100;
    idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idCoupling = -1;
%     idLand = -1; idOS = -1; idRiv = -1; idRad = -1;
    
    problemData.gConst = 10;
    xi0Cont = @(x) zeros(size(x));
    dxZb = 0.005;
    omega = 0.1;
    delta = 0.1;
    epsilon = 0.01;
    rho = 0.001;
    
    xiCont = @(t,x) epsilon * sin(omega * (x+t));
    zBotCont = @(x) -2 + dxZb * x;
    
    problemData.hCont = @(t,x) xiCont(t,x) - zBotCont(x);
    problemData.u1Cont = @(t,x,z) delta * (z - zBotCont(x)) .* sin(omega * (x+t));
    problemData.u2Cont = @(t,x,z) delta * dxZb * z .* sin(omega * (x+t)) - 0.5 * delta * omega * (z - zBotCont(x)).^2 .* cos(omega * (x+t));
    
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
                        
    problemData.DCont = { @(t,x,z) rho * ones(size(x)), @(t,x,z) zeros(size(x)); ...
                          @(t,x,z) zeros(size(x)), @(t,x,z) rho * ones(size(x)) };
                        
    dxzDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
                 @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
    
    problemData.fhCont = @(t,x) dtHCont(t,x) + dxU1hCont(t,x);
    problemData.fuCont = @(t,x,z) dtU1Cont(t,x,z) + 2 * problemData.u1Cont(t,x,z) .* dxU1Cont(t,x,z) + ...
                            dzU1Cont(t,x,z) .* problemData.u2Cont(t,x,z) + problemData.u1Cont(t,x,z) .* dzU2Cont(t,x,z) - ( ...
                              problemData.DCont{1,1}(t,x,z) .* dxdxU1Cont(t,x,z) + dxzDCont{1,1}(t,x,z) .* dxU1Cont(t,x,z) + ...
                              problemData.DCont{1,2}(t,x,z) .* dxdzU1Cont(t,x,z) + dxzDCont{1,2}(t,x,z) .* dzU1Cont(t,x,z) + ...
                              problemData.DCont{2,1}(t,x,z) .* dxdzU1Cont(t,x,z) + dxzDCont{2,1}(t,x,z) .* dxU1Cont(t,x,z) + ...
                              problemData.DCont{2,2}(t,x,z) .* dzdzU1Cont(t,x,z) + dxzDCont{2,2}(t,x,z) .* dzU1Cont(t,x,z) ) + ...
                              problemData.gConst * dxXiCont(t,x);
    
    problemData.hDCont = problemData.hCont;
    problemData.u1DCont = problemData.u1Cont;
    problemData.u2DCont = problemData.u2Cont;
    problemData.q1DCont = @(t,x,z) -problemData.DCont{1,1}(t,x,z) .* dxU1Cont(t,x,z) - problemData.DCont{1,2}(t,x,z) .* dzU1Cont(t,x,z);
    problemData.q2DCont = @(t,x,z) -problemData.DCont{2,1}(t,x,z) .* dxU1Cont(t,x,z) - problemData.DCont{2,2}(t,x,z) .* dzU1Cont(t,x,z);
    problemData.uhDCont = @(t,x) u1hCont(t,x);
  
  case 'coupling'
    domainWidth = 100;
    idLand = [2,4]; idOS = [2,4]; idRiv = [2,4]; idRad = -1; idCoupling = 1;
%     idLand = -1; idOS = -1; idRiv = -1; idRad = -1;
    
    problemData.gConst = 10;
    xi0Cont = @(x) 3 * ones(size(x));
    dxZb = 0;
    a = 1;
    b = 0.01;
    c = 1;
    d = 0.01;
    f = 1;
    k = 1;
    
    xiCont = @(t,x) 3 + a * cos(b*x+c*t);
    zBotCont = @(x) 0 + dxZb * x;
    
    problemData.hCont = @(t,x) xiCont(t,x) - zBotCont(x);
    problemData.u1Cont = @(t,x,z) -k * a * d^2 / b * sin(b*x+c*t) .* cos(d*z) + f * z;
    problemData.u2Cont = @(t,x,z) k * a * d * cos(b*x+c*t) .* sin(d*z);
    
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
                        
    problemData.DCont = cellfun(@(c) @(t,x,z) c * ones(size(x)), {0.01, 0; 0, 0.01}, 'UniformOutput', false);
                        
    dxzDCont = cellfun(@(c) @(t,x,z) c * ones(size(x)), {0, 0; 0, 0}, 'UniformOutput', false);
    
    problemData.fhCont = @(t,x) dtHCont(t,x) + dxU1hCont(t,x);
    problemData.fuCont = @(t,x,z) dtU1Cont(t,x,z) + 2 * problemData.u1Cont(t,x,z) .* dxU1Cont(t,x,z) + ...
                            dzU1Cont(t,x,z) .* problemData.u2Cont(t,x,z) + problemData.u1Cont(t,x,z) .* dzU2Cont(t,x,z) - ( ...
                              problemData.DCont{1,1}(t,x,z) .* dxdxU1Cont(t,x,z) + dxzDCont{1,1}(t,x,z) .* dxU1Cont(t,x,z) + ...
                              problemData.DCont{1,2}(t,x,z) .* dxdzU1Cont(t,x,z) + dxzDCont{1,2}(t,x,z) .* dzU1Cont(t,x,z) + ...
                              problemData.DCont{2,1}(t,x,z) .* dxdzU1Cont(t,x,z) + dxzDCont{2,1}(t,x,z) .* dxU1Cont(t,x,z) + ...
                              problemData.DCont{2,2}(t,x,z) .* dzdzU1Cont(t,x,z) + dxzDCont{2,2}(t,x,z) .* dzU1Cont(t,x,z) ) + ...
                              problemData.gConst * dxXiCont(t,x);
    
    problemData.hDCont = problemData.hCont;
    problemData.u1DCont = problemData.u1Cont;
    problemData.u2DCont = problemData.u2Cont;
    problemData.q1DCont = @(t,x,z) -problemData.DCont{1,1}(t,x,z) .* dxU1Cont(t,x,z) - problemData.DCont{1,2}(t,x,z) .* dzU1Cont(t,x,z);
    problemData.q2DCont = @(t,x,z) -problemData.DCont{2,1}(t,x,z) .* dxU1Cont(t,x,z) - problemData.DCont{2,2}(t,x,z) .* dzU1Cont(t,x,z);
    problemData.uhDCont = @(t,x) u1hCont(t,x);
end % switch
end % function
