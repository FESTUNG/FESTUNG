function hatRoffdiag = computeHatRoffdiagInf(NInf, pd)
p = pd.p;
beta = pd.beta;
gPhiInf1D = pd.basesOnQuadInf.gPhiInf1D;
gThetaPhiInf1D = pd.basesOnQuadInf.gThetaPhiInf1D;
pInf = NInf/(p+1)-1;  qOrdInf = 2*pInf+1; qOrd = 2*p+1;  
hatRoffdiag = zeros(NInf,NInf,NInf,3,3);
[~, W] = quadRule1D(qOrd);
[QInf, WInf] = quadRuleInf1D(qOrdInf, beta);
W = [W, zeros(1, length(WInf)-length(W))];
for nn = 1 : 3
  for np = 1 : 3
    for l = 1 : NInf
      for i = 1 : NInf
        for j = 1 : NInf
          if nn==2 || np==2
            hatRoffdiag(i, j, l, nn,np) = sum( W'    .* gPhiInf1D{qOrdInf}(:,i,nn) .* gThetaPhiInf1D{qOrdInf}(:,l,nn,np) .* gThetaPhiInf1D{qOrdInf}(:,j,nn,np) );
          else
            hatRoffdiag(i, j, l, nn,np) = sum( WInf' .* gPhiInf1D{qOrdInf}(:,i,nn) .* gThetaPhiInf1D{qOrdInf}(:,l,nn,np) .* gThetaPhiInf1D{qOrdInf}(:,j,nn,np) .* exp(-beta*QInf)');
          end % if
        end % for
      end % for
    end % for
  end % for
end % for
end % function
