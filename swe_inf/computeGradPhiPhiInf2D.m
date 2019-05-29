function gradPhiPhi2D = computeGradPhiPhiInf2D(NInf, pd)
p = pd.p;
beta = pd.beta;
gPhiInf2D = pd.basesOnQuadInf.gPhiInf2D;
gGradPhiInf2D = pd.basesOnQuadInf.gGradPhiInf2D;
qOrd = max(2*p, 1);
pInf = NInf/(p+1)-1; qOrdInf = max(2*pInf, 1); [Q1, ~, W] = quadRuleInf2D(qOrd, qOrdInf, beta);
gradPhiPhi2D = zeros(NInf, NInf, length(W), 2);
for i = 1 : NInf
  for j = 1 : NInf
    gradPhiPhi2D(i, j, :, 1) = gGradPhiInf2D{qOrdInf}(:, i, 1) .* gPhiInf2D{qOrdInf}(:, j) .* W.' .* exp(-beta*Q1)';
    gradPhiPhi2D(i, j, :, 2) = gGradPhiInf2D{qOrdInf}(:, i, 2) .* gPhiInf2D{qOrdInf}(:, j) .* W.' .* exp(-beta*Q1)';
  end % for
end % for
end % function
