function hatH = computeHatHInf(NInf, pd)
p = pd.p;
beta = pd.beta;
gPhiInf2D = pd.basesOnQuadInf.gPhiInf2D;
gGradPhiInf2D = pd.basesOnQuadInf.gGradPhiInf2D;
qOrd = max(2*p, 1);
pInf = NInf/(p+1)-1; qOrdInf = max(2*pInf, 1); [Q1, ~, W] = quadRuleInf2D(qOrd, qOrdInf, beta);
hatH = zeros(NInf, NInf, 2, length(W));
for i = 1 : NInf
  for j = 1 : NInf
    for m = 1 : 2
      hatH(i, j, m, :) = W.' .* gGradPhiInf2D{qOrdInf}(:,i,m) .* gPhiInf2D{qOrdInf}(:,j) .* exp(-beta*Q1)';
    end % for
  end % for
end % for
end % function
