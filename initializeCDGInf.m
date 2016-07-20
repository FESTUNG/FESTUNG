function cDGInf = initializeCDGInf(gInf, xi0, u0, v0, zbInfQuad, ord, ordInf, hatM, beta, basesOnQuad)
% global gPhiInf2D
ordInf = max(ordInf,1);  [Q1, Q2, W] = quadRuleInf2D(ord, ordInf, beta);
NInf = size(hatM, 1);
F1 = @(X1, X2) gInf.xHat(:,1)*ones(size(X1)) + cos(gInf.Alpha)*X1 - ...
  (gInf.H*ones(size(X1)) + 2*gInf.M*X1).* (sin(gInf.Beta)*(X2-0.5));
F2 = @(X1, X2) gInf.xHat(:,2)*ones(size(X2)) + sin(gInf.Alpha)*X1 + ...
  (gInf.H*ones(size(X1)) + 2*gInf.M*X1).* (cos(gInf.Beta)*(X2-0.5));
%note: det is NOT constant => projection element-wise 
rhs = zeros(NInf,1);
cDGInf = zeros(gInf.numInfElem, NInf, 3);
for l = 1 : 3
  switch l
    case 1
      c = xi0(F1(Q1, Q2), F2(Q1, Q2))-zbInfQuad;
    case 2
      c = (xi0(F1(Q1, Q2), F2(Q1, Q2))-zbInfQuad) .* u0(F1(Q1, Q2), F2(Q1, Q2));
    case 3
      c = (xi0(F1(Q1, Q2), F2(Q1, Q2))-zbInfQuad) .* v0(F1(Q1, Q2), F2(Q1, Q2));
  end %switch
  for k = 1 : gInf.numInfElem
    hatMk = zeros(NInf);
    for i = 1 : NInf
      rhs(i) = sum (c(k,:).*gInf.detQuadPoints(k,:).*basesOnQuad.gPhiInf2D{ordInf}(:,i)'.*W.*exp(-beta*Q1));
    end % for
    for r = 1 : size(Q1,2)
      hatMk = hatMk + hatM(:,:,r)*gInf.detQuadPoints(k,r);
    end %for
    cDGInf(k,:,l) = (hatMk \ rhs)';
  end % for
end % for
end % function
