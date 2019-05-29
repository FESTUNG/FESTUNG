function globH = assembleGlobHInf(gInf, hatH, p, beta)
%% assembles global contributions of - \int_T (\nabla \varphi_i) * \varphi_j
%% multiplying globH with uH, vH corresponds to area integrals in weak form of div( [uH, vH]^T ) in PDE

NInf = size(hatH,2); KInf = gInf.numInfElem;
qOrd = max(2*p, 1);
pInf = NInf/(p+1)-1; qOrdInf = max(2*pInf, 1); [Q1,Q2,~] = quadRuleInf2D(qOrd, qOrdInf, beta);
globH = cell(2, 1);
globH{1} = sparse(KInf*NInf, KInf*NInf);  globH{2} = sparse(KInf*NInf, KInf*NInf);
for r = 1 : size(Q1,2)
  globH{1} = globH{1} - kron(spdiags(cos(gInf.Beta).*(gInf.H+2*Q1(r)*gInf.M), 0, KInf, KInf), hatH(:,:,1,r)) ...
                      - kron(spdiags(-sin(gInf.Alpha)-cos(gInf.Beta)*2.*gInf.M*(Q2(r)-0.5), 0, KInf, KInf), hatH(:,:,2,r));
  globH{2} = globH{2} - kron(spdiags(sin(gInf.Beta).*(gInf.H+2*Q1(r)*gInf.M), 0, KInf, KInf), hatH(:,:,1,r)) ...
                      - kron(spdiags( cos(gInf.Alpha)-sin(gInf.Beta)*2.*gInf.M*(Q2(r)-0.5), 0, KInf, KInf), hatH(:,:,2,r));
end%for
end % function
