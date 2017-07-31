function problemData = getTestcase(problemData, testcase)  
% checkMultipleIds = @(idE0T, ids) logical(sum(bsxfun(@eq, idE0T, reshape(ids, 1, 1, length(ids))), 3));
switch testcase
  case 'solid_body' % LeVeque's solid body rotation
    isAnalytical = false;
    isStationary = false;
    
    G = @(x1, x2, x1_0, x2_0) (1/0.15) * sqrt((x1-x1_0).^2 + (x2-x2_0).^2);
    c0Cont = @(x1, x2) ((x1 - 0.5).^2 + (x2 - 0.75).^2 <= 0.0225 & (x1 <= 0.475 | x1 >= 0.525 | x2 >= 0.85)) + ...
        (1-G(x1, x2, 0.5, 0.25)) .* ((x1 - 0.5).^2 + (x2 - 0.25).^2 <= 0.0225) + ...
        0.25*(1+cos(pi*G(x1, x2, 0.25, 0.5))).*((x1 - 0.25).^2 + (x2 - 0.5).^2 <= 0.0225);
    fCont = @(t,x1,x2) zeros(size(x1));
    u1Cont = @(t,x1,x2) 0.5 - x2;
    u2Cont = @(t,x1,x2) x1 - 0.5;
    cDCont = @(t,x1,x2) zeros(size(x1));

    numSteps = 320;
    tEnd = 2 * pi;

  case 'stationary' % Stationary analytical example
    isAnalytical = true;
    isStationary = true;
    
    cCont = @(t, x1, x2) cos(7 * x1) .* cos(7 * x2);
    u1Cont = @(t, x1, x2) exp(0.5 * (x1 + x2));
    u2Cont = @(t, x1, x2) exp(0.5 * (x1 - x2));
    fCont = @(t, x1, x2) -7 * u1Cont(t, x1, x2) .* sin(7 * x1) .* cos(7 * x2) ...
        - 7 * u2Cont(t, x1, x2) .* cos(7 * x1) .* sin(7 * x2) ...
        + 0.5 * (u1Cont(t, x1, x2) - u2Cont(t, x1, x2)) .* cCont(t, x1, x2);
    c0Cont = @(x1, x2) cCont(0, x1, x2);
    cDCont = @(t, x1, x2) cCont(t, x1, x2);
    
    numSteps = 1;
    tEnd = 10;
    
  case 'transient_ode' % Transient ODE example
    isAnalytical = true;
    isStationary = false;
    
    cCont = @(t, x1, x2) exp(-t) * ones(size(x1));
    u1Cont = @(t, x1, x2) zeros(size(x1));
    u2Cont = @(t, x1, x2) zeros(size(x1));
    fCont = @(t, x1, x2) -cCont(t, x1, x2);
    c0Cont = @(x1, x2) cCont(0, x1, x2);
    cDCont = @(t, x1, x2) cCont(t, x1, x2);
    
    numSteps = 4;
    tEnd = 2;
    
  case 'transient' % Transient analytical example
    isAnalytical = true;
    isStationary = false;
    
    cCont = @(t, x1, x2) cos(7 * x1) .* cos(7 * x2) + exp(-t);
    u1Cont = @(t, x1, x2) exp(0.5 * (x1 + x2));
    u2Cont = @(t, x1, x2) exp(0.5 * (x1 - x2));
    fCont = @(t, x1, x2) -exp(-t) - 7 * u1Cont(t, x1, x2) .* sin(7 * x1) .* cos(7 * x2) ...
        - 7 * u2Cont(t, x1, x2) .* cos(7 * x1) .* sin(7 * x2) ...
        + 0.5 * (u1Cont(t, x1, x2) - u2Cont(t, x1, x2)) .* cCont(t, x1, x2);
    c0Cont = @(x1, x2) cCont(0, x1, x2);
    cDCont = @(t, x1, x2) cCont(t, x1, x2);
    
    numSteps = 10;
    tEnd = 2;
   
  otherwise
    error('Invalid testcase "%s".', testcase);
end % switch

fprintf('Loaded testcase "%s".\n', testcase);

% Time stepping parameters
problemData = setdefault(problemData, 'numSteps', numSteps);
problemData = setdefault(problemData, 'tEnd', tEnd);
problemData = setdefault(problemData, 'isStationary', isStationary);

% Store defined functions in struct
problemData = setdefault(problemData, 'c0Cont', c0Cont);
problemData = setdefault(problemData, 'fCont', fCont);
problemData = setdefault(problemData, 'u1Cont', u1Cont);
problemData = setdefault(problemData, 'u2Cont', u2Cont);
problemData = setdefault(problemData, 'cDCont', cDCont);
if isAnalytical, problemData = setdefault(problemData, 'cCont', cCont); end

% Specify boundary conditions
problemData = setdefault(problemData, 'generateMarkE0Tint', @(g) g.idE0T == 0);
problemData = setdefault(problemData, 'generateMarkE0Tbdr', @(g) g.idE0T ~= 0);

end % function
