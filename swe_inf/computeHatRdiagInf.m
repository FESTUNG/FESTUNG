function hatRdiag = computeHatRdiagInf(NInf, pd)
p = pd.p;
beta = pd.beta;
gPhiInf1D = pd.basesOnQuadInf.gPhiInf1D;
pInf = NInf/(p+1)-1; qOrdInf = 2*pInf+1; qOrd = 2*p+1;
hatRdiag = zeros(NInf, NInf, NInf,3);
[~, W] = quadRule1D(qOrd);
[QInf, WInf] = quadRuleInf1D(qOrdInf, beta);
W = [W, zeros(1, length(WInf)-length(W))];
for n = 1 : 3
  for l = 1 : NInf
    for i = 1 : NInf
      for j = 1 : NInf
        if n == 2
          hatRdiag(i,j,l,n) = sum(W'    .* gPhiInf1D{qOrdInf}(:,i,n) .* gPhiInf1D{qOrdInf}(:,l,n) .* gPhiInf1D{qOrdInf}(:,j,n));
        else
          hatRdiag(i,j,l,n) = sum(WInf' .* gPhiInf1D{qOrdInf}(:,i,n) .* gPhiInf1D{qOrdInf}(:,l,n) .* gPhiInf1D{qOrdInf}(:,j,n).*exp(-beta*QInf)');          
        end%if
      end % for
    end % for
  end % for
end % for
end % function
