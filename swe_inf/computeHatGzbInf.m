%ohne aufsummieren!
function hatGzb = computeHatGzbInf(NInf, pd)
p = pd.p;
beta = pd.beta;
gPhiInf2D = pd.basesOnQuadInf.gPhiInf2D;
qOrd = max(2*p, 1);
pInf = NInf/(p+1)-1; qOrdInf = max(2*pInf, 1); [Q1, Q2, W] = quadRuleInf2D(qOrd, qOrdInf, beta);
hatGzb = zeros(NInf, NInf, 2, 2, size(Q1,2));
for i = 1 : NInf
  for j = 1 : NInf
    for l = 1 : 2%2 lin Ansatzfunktionen
      for m = 1 : 2
        hatGzb(i,j,l,m,:) = W.' .* gradPhiLin(l, m, Q1, Q2).' .* gPhiInf2D{qOrdInf}(:,i) .* gPhiInf2D{qOrdInf}(:,j) .* exp(-beta*Q1)';
      end % for
    end % for
  end % for
end % for
end % function