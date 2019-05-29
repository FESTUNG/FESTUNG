function globFExp = assembleGlobFExpInf(gInf, hatD, cDG, gConst, p, beta)
%% assembles global contributions of 0.5 * beta * gConst * \int_T (\varphi_i) * H * \varphi_j
%% multiplying globF with H corresponds to area integrals in weak form of 0.5 * gConst * H^2 in PDE
[KInf, NInf, ~] = size(cDG);
qOrd = max(2*p,1);
pInf = NInf/(p+1)-1; qOrdInf = max(2*pInf,1);
[Q1,~,~] = quadRuleInf2D(qOrd, qOrdInf, beta);
globFExp = cell(2, 1);
globFExp{1} = sparse(KInf*NInf, KInf*NInf); globFExp{2} = sparse(KInf*NInf, KInf*NInf);
for l = 1 : NInf
  for r = 1 : size(Q1,2)
    globFExp{1} = globFExp{1} + 0.5* beta * gConst * kron(spdiags(cDG(:,l,1).*cos(gInf.Beta).*(gInf.H+2*Q1(r)*gInf.M), 0, KInf, KInf), hatD(:,:,l,r));
    globFExp{2} = globFExp{2} + 0.5* beta * gConst * kron(spdiags(cDG(:,l,1).*sin(gInf.Beta).*(gInf.H+2*Q1(r)*gInf.M), 0, KInf, KInf), hatD(:,:,l,r));
  end % for
end % for
end % function