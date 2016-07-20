function globHExp = assembleGlobHExpInf(gInf, hatM, p, beta)
%% assembles global contributions of beta * \int_T \varphi_i * \varphi_j

NInf = size(hatM,2); KInf = gInf.numInfElem;
qOrd = max(2*p, 1);
pInf = NInf/(p+1)-1; qOrdInf = max(2*pInf, 1); [Q1,~,~] = quadRuleInf2D(qOrd, qOrdInf, beta);
globHExp = cell(2, 1);
globHExp{1} = sparse(KInf*NInf, KInf*NInf);  globHExp{2} = sparse(KInf*NInf, KInf*NInf);
for r = 1 : size(Q1,2)
  globHExp{1} = globHExp{1} + beta * kron(spdiags(cos(gInf.Beta).*(gInf.H+2*Q1(r)*gInf.M), 0, KInf, KInf), hatM(:,:,r));
  globHExp{2} = globHExp{2} + beta * kron(spdiags(sin(gInf.Beta).*(gInf.H+2*Q1(r)*gInf.M), 0, KInf, KInf), hatM(:,:,r));                            
end%for
end % function
