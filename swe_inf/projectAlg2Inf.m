function cInf = projectAlg2Inf(gInf, cAlg, ord, ordInf, hatM, beta, basesOnQuad)
% global gPhiInf2D
ord = max(ord,1); ordInf = max(ordInf,1);  [Q1, Q2, W] = quadRuleInf2D(ord, ordInf, beta);
NInf = size(hatM, 1);

F1 = @(X1, X2) gInf.xHat(:,1)*ones(size(X1)) + cos(gInf.Alpha)*X1 - ...
  (gInf.H*ones(size(X1)) + 2*gInf.M*X1).* (sin(gInf.Beta)*(X2-0.5));
F2 = @(X1, X2) gInf.xHat(:,2)*ones(size(X2)) + sin(gInf.Alpha)*X1 + ...
  (gInf.H*ones(size(X1)) + 2*gInf.M*X1).* (cos(gInf.Beta)*(X2-0.5));

rhs = zeros(NInf,1);
cInf = zeros(gInf.numInfElem, NInf);
c = cAlg(F1(Q1, Q2), F2(Q1, Q2));
for k = 1 : gInf.numInfElem
  hatMk = zeros(NInf);
  for i = 1 : NInf
    rhs(i) = sum (c(k,:).*gInf.detQuadPoints(k,:).*basesOnQuad.gPhiInf2D{ordInf}(:,i)'.*W.*exp(-beta*Q1));
  end % for
  for r = 1 : size(Q1,2)
    hatMk = hatMk + hatM(:,:,r).*gInf.detQuadPoints(k,r);
  end %for
  cInf(k,:) = (hatMk \ rhs)';
end % for
end % function
