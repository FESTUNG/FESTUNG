function PhiPhi2D = computePhiPhiInf2D(NInf, pd)
p = pd.p;
beta = pd.beta;
gPhiInf2D = pd.basesOnQuadInf.gPhiInf2D;
qOrd = max(2*p, 1);
pInf = NInf/(p+1)-1; qOrdInf = max(2*pInf, 1); [Q1, ~, W] = quadRuleInf2D(qOrd, qOrdInf, beta);
PhiPhi2D = zeros(NInf, NInf, length(W));
for i = 1 : NInf
  for j = 1 : NInf
    PhiPhi2D(i, j, :) = gPhiInf2D{qOrdInf}(:, i) .* gPhiInf2D{qOrdInf}(:, j) .* W.' .* exp(-beta*Q1)';
  end % for
end % for
end % function
