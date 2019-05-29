function UDG = getInfvelocities(gInf, cDG, ord, ordInf, hatM, beta)
%% approximates the velocities by projection of quotients of primary unknowns
global gPhiInf2D
ordInf = max(ordInf,1);  [Q1, ~, W] = quadRuleInf2D(ord, ordInf, beta);
KInf = gInf.numInfElem; NInf = size(hatM, 1);
UDG = zeros(KInf,NInf,2);

rhs1 = zeros(NInf,1);
rhs2 = zeros(NInf,1);
c1 = ( cDG(:,:,2) * gPhiInf2D{ordInf}(:,:).' ) ./ ( cDG(:,:,1) * gPhiInf2D{ordInf}(:,:).' );
c2 = ( cDG(:,:,3) * gPhiInf2D{ordInf}(:,:).' ) ./ ( cDG(:,:,1) * gPhiInf2D{ordInf}(:,:).' );
for k = 1 : gInf.numInfElem
  hatMk = zeros(NInf);
  for i = 1 : NInf
    rhs1(i) = sum (c1(k,:).*gInf.detQuadPoints(k,:).*gPhiInf2D{ordInf}(:,i)'.*W.*exp(-beta*Q1));
    rhs2(i) = sum (c2(k,:).*gInf.detQuadPoints(k,:).*gPhiInf2D{ordInf}(:,i)'.*W.*exp(-beta*Q1));
  end % for
  for r = 1 : size(W,2)
    hatMk = hatMk + hatM(:,:,r)*gInf.detQuadPoints(k,r);
  end %for
  UDG(k,:,1) = (hatMk \ rhs1)';
  UDG(k,:,2) = (hatMk \ rhs2)';
end % for
end % function
