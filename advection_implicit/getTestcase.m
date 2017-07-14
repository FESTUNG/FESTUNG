function problemData = getTestcase(problemData, testcase)
checkMultipleIds = @(idE0T, ids) logical(sum(bsxfun(@eq, idE0T, reshape(ids, 1, 1, length(ids))), 3));
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
    gNCont = @(t,x1,x2) zeros(size(x1));
    
    generateMarkE0TbdrN = @(g) false(g.numT, 3);
    
    numSteps = 160;
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
    gNCont = @(t, x1, x2) -7 * cos(7 * x1) .* sin(7 * x2);
    
    idBdrN = -1;
    generateMarkE0TbdrN = @(g) checkMultipleIds(g.idE0T, idBdrN);
    
    numSteps = 1;
    tEnd = 10;
    
  case 'transient' % Transient analytical example
    isAnalytical = true;
    isStationary = false;
    
    cCont = @(t, x1, x2) cos(7 * x1) .* cos(7 * x2) + t;
    u1Cont = @(t, x1, x2) exp(0.5 * (x1 + x2));
    u2Cont = @(t, x1, x2) exp(0.5 * (x1 - x2));
    fCont = @(t, x1, x2) 1 - 7 * u1Cont(t, x1, x2) .* sin(7 * x1) .* cos(7 * x2) ...
        - 7 * u2Cont(t, x1, x2) .* cos(7 * x1) .* sin(7 * x2) ...
        + 0.5 * (u1Cont(t, x1, x2) - u2Cont(t, x1, x2)) .* cCont(t, x1, x2);
    c0Cont = @(x1, x2) cCont(0, x1, x2);
    cDCont = @(t, x1, x2) cCont(t, x1, x2);
    gNCont = @(t, x1, x2) zeros(size(x1));
    
    generateMarkE0TbdrN = @(g) false(g.numT, 3);
                              
    numSteps = 10;
    tEnd = 20;    
    
  case 'transient2' % Transient analytical example
    isAnalytical = true;
    isStationary = false;
    
    cCont = @(t, x1, x2) cos(7 * x1) .* cos(7 * x2) + sin(t);
    u1Cont = @(t, x1, x2) exp(0.5 * (x1 + x2));
    u2Cont = @(t, x1, x2) exp(0.5 * (x1 - x2));
    fCont = @(t, x1, x2) cos(t) - 7 * u1Cont(t, x1, x2) .* sin(7 * x1) .* cos(7 * x2) ...
        - 7 * u2Cont(t, x1, x2) .* cos(7 * x1) .* sin(7 * x2) ...
        + 0.5 * (u1Cont(t, x1, x2) - u2Cont(t, x1, x2)) .* cCont(t, x1, x2);
    c0Cont = @(x1, x2) cCont(0, x1, x2);
    cDCont = @(t, x1, x2) cCont(t, x1, x2);
    gNCont = @(t, x1, x2) -7 * cos(7 * x1) .* sin(7 * x2);
    
    generateMarkE0TbdrN = @(g) false(g.numT, 3);
                              
    numSteps = 10;
    tEnd = 40;    
   
  otherwise
    error('Invalid testcase "%s".', testcase);
end % switch

fprintf('Loaded testcase "%s".\n', testcase);

% Time stepping parameters
problemData = setdefault(problemData, 'isStationary', isStationary);
problemData = setdefault(problemData, 'numSteps', numSteps);
problemData = setdefault(problemData, 'tEnd', tEnd);

% Store defined functions in struct
problemData = setdefault(problemData, 'c0Cont', c0Cont);
problemData = setdefault(problemData, 'fCont', fCont);
problemData = setdefault(problemData, 'u1Cont', u1Cont);
problemData = setdefault(problemData, 'u2Cont', u2Cont);
problemData = setdefault(problemData, 'cDCont', cDCont);
problemData = setdefault(problemData, 'gNCont', gNCont);
if isAnalytical, problemData = setdefault(problemData, 'cCont', cCont); end

% Specify boundary conditions
problemData = setdefault(problemData, 'generateMarkE0Tint', @(g) g.idE0T == 0);
problemData = setdefault(problemData, 'generateMarkE0TbdrN', generateMarkE0TbdrN);
problemData = setdefault(problemData, 'generateMarkE0TbdrD', @(g) ~(g.markE0Tint | g.markE0TbdrN));
end % function