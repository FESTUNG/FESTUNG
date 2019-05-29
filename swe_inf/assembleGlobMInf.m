function globM = assembleGlobMInf(gInf, hatM, p, beta)
%% assembles global contributions of \int_T \varphi_i * \varphi_j
%% globM is the mass matrix
NInf = size(hatM,2); KInf = gInf.numInfElem;
qOrd = max(2*p, 1);
pInf = NInf/(p+1)-1; qOrdInf = max(2*pInf, 1); [Q1,~,~] = quadRuleInf2D(qOrd, qOrdInf, beta);
globM = sparse(KInf*NInf,KInf*NInf);
for r = 1 : size(Q1,2)
  globM = globM + kron(spdiags(gInf.detQuadPoints(:,r), 0, KInf, KInf), hatM(:,:,r));
end % for
end % function
