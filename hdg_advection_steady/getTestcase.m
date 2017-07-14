function problemData = getTestcase(problemData, problemName)
switch problemName
  case 'Convergence'
    isAnalytical = true;
    cCont  = @(x1,x2) cos(7*x1).*cos(7*x2);
    c0Cont = @(x1,x2) cCont(x1,x2);
    cDCont = @(x1,x2) cCont(x1,x2);
    fCont  = @(x1,x2) -7*sin(7*x1).*cos(7*x2).*exp((x1+x2)/2) ...
      -7*cos(7*x1).*sin(7*x2).*exp((x1-x2)/2) ...
      + 0.5*cos(7*x1).*cos(7*x2).*exp((x1+x2)/2) ...
      - 0.5*cos(7*x1).*cos(7*x2).*exp((x1-x2)/2);
%     fluxCont = @(x1, x2, c ) evalSteadyFlux(0, x1, x2, c);
    u1Cont = @(x1,x2) exp((x1+x2)/2);
    u2Cont = @(x1,x2) exp((x1-x2)/2);
    
%     generateMarkE0TbdrN = @(g) generateSteadyOutflowBoundary(g);
%     generateMarkE0TbdrD = @(g) ~(g.markE0Tint | g.markE0TbdrN);
    idBdrN = [2 3];
    
    generateGridData = @(hmax) domainSquare(hmax);
    
  otherwise
    error('Specified test case is not available. Please check your configuration.');


end % switch

problemData = setdefault(problemData, 'c0Cont', c0Cont);
problemData = setdefault(problemData, 'u1Cont', u1Cont);
problemData = setdefault(problemData, 'u2Cont', u2Cont);
problemData = setdefault(problemData, 'cCont', cCont);
problemData = setdefault(problemData, 'fCont', fCont);
% problemData = setdefault(problemData, 'fluxCont', fluxCont);
problemData = setdefault(problemData, 'cDCont', cDCont);
if isAnalytical, problemData = setdefault(problemData, 'cCont', cCont); end

% Specify triangulation and edge ids of boundary conditions
problemData = setdefault(problemData, 'generateGridData', generateGridData);

checkMultipleIds = @(idE0T, ids) logical(sum(bsxfun(@eq, idE0T, reshape(ids, 1, 1, length(ids))), 3));
problemData = setdefault(problemData, 'generateMarkE0Tint', @(g) g.idE0T == 0);
problemData = setdefault(problemData, 'generateMarkE0TbdrN', @(g) checkMultipleIds(g.idE0T, idBdrN));
problemData = setdefault(problemData, 'generateMarkE0TbdrD', @(g) ~(g.markE0Tint | g.markE0TbdrN));

end % function
