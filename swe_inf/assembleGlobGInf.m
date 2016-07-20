function globG = assembleGlobGInf(gInf, hatGzb, zbLin, gConst, p, beta)%hatGzb noch NICHT aufsummiert!
%% assembles global contributions of gConst * \int_T \varphi_i * (\nabla z_b) * \varphi_j
%% multiplying globG with H corresponds to area integrals in weak form of gConst * (\nabla z_b) * H in PDE
NInf = size(hatGzb,2); KInf = size(zbLin,1);
qOrd = max(2*p, 1);
pInf = NInf/(p+1)-1; qOrdInf = max(2*pInf, 1); [Q1,Q2,~] = quadRuleInf2D(qOrd, qOrdInf, beta);
globG = cell(2, 1);
globG{1} = sparse(KInf*NInf, KInf*NInf);  globG{2} = sparse(KInf*NInf, KInf*NInf);
for l = 1 : 2%geaendert
  for r = 1 : size(Q1,2)
    globG{1} = globG{1} + gConst * kron(spdiags(zbLin(:,l).*  cos(gInf.Beta).*(gInf.H+2*Q1(r)*gInf.M)              , 0, KInf, KInf), hatGzb(:,:,l,1,r)) ...
                        + gConst * kron(spdiags(zbLin(:,l).*(-sin(gInf.Alpha)-cos(gInf.Beta)*2.*gInf.M*(Q2(r)-0.5)), 0, KInf, KInf), hatGzb(:,:,l,2,r));
    globG{2} = globG{2} + gConst * kron(spdiags(zbLin(:,l).*  sin(gInf.Beta).*(gInf.H+2*Q1(r)*gInf.M)              , 0, KInf, KInf), hatGzb(:,:,l,1,r)) ...
                        + gConst * kron(spdiags(zbLin(:,l).*( cos(gInf.Alpha)-sin(gInf.Beta)*2.*gInf.M*(Q2(r)-0.5)), 0, KInf, KInf), hatGzb(:,:,l,2,r));
  end%for
end % for
end % function
