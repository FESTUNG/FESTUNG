function globF = assembleGlobFInf(gInf, hatG, cDGInf, gConst, p, beta)
%% assembles global contributions of -0.5 * gConst * \int_T (\nabla \varphi_i) * H * \varphi_j
%% multiplying globF with H corresponds to area integrals in weak form of 0.5 * gConst * H^2 in PDE
[KInf, NInf, ~] = size(cDGInf);
qOrd = max(2*p,1);
pInf = NInf/(p+1)-1; qOrdInf = max(2*pInf,1);
[Q1,Q2,~] = quadRuleInf2D(qOrd, qOrdInf, beta);
globF = cell(2, 1);
globF{1} = sparse(KInf*NInf, KInf*NInf); globF{2} = sparse(KInf*NInf, KInf*NInf);
for l = 1 : NInf
  for r = 1 : size(Q1,2)
    globF{1} = globF{1} - 0.5 * gConst * ( kron(spdiags(cDGInf(:,l,1).*cos(gInf.Beta).*(gInf.H+2*Q1(r)*gInf.M), 0, KInf, KInf), hatG(:,:,l,1,r)) ...
                                         + kron(spdiags(cDGInf(:,l,1).*(-sin(gInf.Alpha)-cos(gInf.Beta)*2.*gInf.M*(Q2(r)-0.5)), 0, KInf, KInf), hatG(:,:,l,2,r)) );
    globF{2} = globF{2} - 0.5 * gConst * ( kron(spdiags(cDGInf(:,l,1).*(cos(gInf.Alpha)-sin(gInf.Beta)*2.*gInf.M*(Q2(r)-0.5)), 0, KInf, KInf), hatG(:,:,l,2,r)) ...
                                         + kron(spdiags(cDGInf(:,l,1).*sin(gInf.Beta).*(gInf.H+2*Q1(r)*gInf.M), 0, KInf, KInf), hatG(:,:,l,1,r)) );
  end % for
end % for
end % function