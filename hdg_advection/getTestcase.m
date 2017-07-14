function problemData = getTestcase(problemData, problemName)
switch problemName
  case 'LeVeque'
    isAnalytical = false;
    
    G = @(x1, x2, x1_0, x2_0) (1/0.15) * sqrt((x1-x1_0).^2 + (x2-x2_0).^2);
    %     Initial conditions
    c0Cont = @(x1, x2) ((x1 - 0.5).^2 + (x2 - 0.75).^2 <= 0.0225 & (x1 <= 0.475 | x1 >= 0.525 | x2 >= 0.85)) + ...
      (1-G(x1, x2, 0.5, 0.25)) .* ((x1 - 0.5).^2 + (x2 - 0.25).^2 <= 0.0225) + ...
      0.25*(1+cos(pi*G(x1, x2, 0.25, 0.5))).*((x1 - 0.25).^2 + (x2 - 0.5).^2 <= 0.0225);
    %     Source term
    fCont = @(t,x1,x2) zeros(size(x1));
    %     Flux function
%     fluxCont = @( t, x1, x2, c ) evalLeVequeFlux(t, x1, x2, c);
    u1Cont = @(t, x1, x2) 0.5 - x2;
    u2Cont = @(t, x1, x2) x1 - 0.5;
    %     Dirichlet boundary data
    cDCont = @(t,x1,x2) zeros(size(x1));
    %     Solution
    
%     idBdrN = 
    generateMarkE0TbdrN = @(g) (u1Cont(0, g.baryE0T(:, :, 1), g.baryE0T(:, :, 2)) .* g.nuE0T(:, :, 1) + ...
                                u2Cont(0, g.baryE0T(:, :, 1), g.baryE0T(:, :, 2)) .* g.nuE0T(:, :, 2)) > 0 & ...
                                g.idE0T ~= 0;
    
    
    generateGridData = @(hmax) domainSquare(hmax, 0, 1);

    tEnd = 2*pi;
    
%   case 'RotatingGaussian'
%     getLinearAdvectionSol = @(t, X1, X2) sin( 2. * pi .* (X1 - 1. * t) ) .* sin( 2. * pi .* (X2 - 1. * t)  );
%     c0Cont = @(x1, x2) getLinearAdvectionSol(0, x1, x2);
%     fCont = @(t,x1,x2) zeros(size(x1));
%     fluxCont = @( t, x1, x2, c ) evalRotatingGaussianFlux(t, x1, x2, c);
%     
%     cDCont = @(t,x1,x2) getLinearAdvectionSol(t, x1, x2);
%     gNCont = @(t,x1,x2) zeros(size(x1));
%     
%     generateMarkE0TbdrN = @(g) generateRotGaussBoundary(g);
%     generateMarkE0TbdrD = @(g) ~(g.markE0Tint | g.markE0TbdrN);
%     
%     generateGridData = @(hmax) domainArbitrarySquare( -0.5, 0.5, hmax );
    
  case 'SteadyProblem'
    isAnalytical = true;
    
    cCont  = @(t,x1,x2) cos(7*x1).*cos(7*x2);
    c0Cont = @(x1,x2) cCont(0,x1,x2);
    cDCont = @(t,x1,x2) cCont(t,x1,x2);
    fCont  = @(t, x1,x2) -7*sin(7*x1).*cos(7*x2).*exp((x1+x2)/2) ...
      -7*cos(7*x1).*sin(7*x2).*exp((x1-x2)/2) ...
      + 0.5*cos(7*x1).*cos(7*x2).*exp((x1+x2)/2) ...
      - 0.5*cos(7*x1).*cos(7*x2).*exp((x1-x2)/2);
%     fluxCont = @(t,  x1, x2, c ) evalSteadyFlux(0, x1, x2, c);
    u1Cont = @(t, x1, x2) exp((x1+x2)/2);
    u2Cont = @(t, x1, x2) exp((x1-x2)/2);
    
    idBdrN = [2 3];
%     generateMarkE0TbdrN = @(g) generateSteadyOutflowBoundary(g);
    checkMultipleIds = @(idE0T, ids) logical(sum(bsxfun(@eq, idE0T, reshape(ids, 1, 1, length(ids))), 3));
    generateMarkE0TbdrN = @(g) checkMultipleIds(g.idE0T, idBdrN);
    
    generateGridData = @(hmax) domainSquare(hmax, 0, 1);
  otherwise
    error('Specified test case is not available. Please check your configuration.');


end % switch

problemData = setdefault(problemData, 'c0Cont', c0Cont);
problemData = setdefault(problemData, 'fCont', fCont);
problemData = setdefault(problemData, 'u1Cont', u1Cont);
problemData = setdefault(problemData, 'u2Cont', u2Cont);
problemData = setdefault(problemData, 'cDCont', cDCont);
if isAnalytical, problemData = setdefault(problemData, 'cCont', cCont); end

% Specify triangulation and edge ids of boundary conditions
problemData = setdefault(problemData, 'generateGridData', generateGridData);

checkMultipleIds = @(idE0T, ids) logical(sum(bsxfun(@eq, idE0T, reshape(ids, 1, 1, length(ids))), 3));
problemData = setdefault(problemData, 'generateMarkE0Tint', @(g) g.idE0T == 0);
problemData = setdefault(problemData, 'generateMarkE0TbdrN', generateMarkE0TbdrN);%@(g) checkMultipleIds(g.idE0T, idbdrN));
problemData = setdefault(problemData, 'generateMarkE0TbdrD', @(g) ~(g.markE0Tint | g.markE0TbdrN));




end % function
