function globD = assembleGlobDInf(gInf, hatD, fcInf, p, beta)
%% assembles global contributions of \int_T \varphi_i * fc * \varphi_j
%% multiplying globD with uH, vH corresponds to area integrals in weak form of fc * uH, fc * vH respectively in PDE
KInf = gInf.numInfElem; NInf= size(hatD,1);
qOrd = max(2*p, 1);
pInf = NInf/(p+1)-1; qOrdInf = max(2*pInf, 1); [Q1,~,~] = quadRuleInf2D(qOrd, qOrdInf, beta);
globD = sparse(KInf*NInf,KInf*NInf);
for l = 1 : NInf
  for r = 1 : size(Q1,2)
    globD = globD + kron(spdiags(gInf.detQuadPoints(:,r) .* fcInf(:,l), 0, KInf, KInf), hatD(:,:,l,r));
  end % for
end % for
end % function
