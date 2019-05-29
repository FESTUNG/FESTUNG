function PhiPhidiag = computePhiPhidiagInf(NInf, pd)
p = pd.p;
beta = pd.beta;
gPhiInf1D = pd.basesOnQuadInf.gPhiInf1D;
pInf = NInf/(p+1)-1;  qOrdInf = max(2*pInf+1, 1); qOrd = 2*p+1;  [~, W] = quadRuleInf1D(qOrdInf, beta);
PhiPhidiag = zeros(NInf, NInf, length(W) ,3);
[~, W] = quadRule1D(qOrd);
[QInf, WInf] = quadRuleInf1D(qOrdInf, beta);
W = [W, zeros(1, length(WInf)-length(W))];
for n = 1 : 3
  for i = 1 : NInf
    for j = 1 : NInf
      if n == 2
        PhiPhidiag(i, j, :, n) = gPhiInf1D{qOrdInf}(:, i, n) .* gPhiInf1D{qOrdInf}(:, j, n) .* W.';
      else
        PhiPhidiag(i, j, :, n) = gPhiInf1D{qOrdInf}(:, i, n) .* gPhiInf1D{qOrdInf}(:, j, n) .* WInf.' .* exp(-beta*QInf)';
      end % if
    end % for
  end % for
end % for
end % function
