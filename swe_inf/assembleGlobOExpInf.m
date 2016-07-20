function globOExp = assembleGlobOExpInf(gInf, PhiPhi2D, cDG, p, beta, gC2DInf)
%% assembles global contributions of beta * \int_T ( ( \varphi_i ) * uH / H ) * \varphi_j
% global gC2DInf
[KInf, NInf, ~] = size(cDG);
qOrd = max(2*p,1);
pInf = NInf/(p+1)-1;  qOrdInf = max(2*pInf, 1); [Q1, ~, ~] = quadRuleInf2D(qOrd, qOrdInf, beta);
globOExp = cell(2, 1);
globOExp{1} = sparse(KInf*NInf, KInf*NInf); globOExp{2} = sparse(KInf*NInf, KInf*NInf);
for r = 1 : length(Q1)
  globOExp{1} = globOExp{1} + beta * kron( spdiags( cos(gInf.Beta).*(gInf.H+2*Q1(r)*gInf.M) .* gC2DInf{2}(:,r)  ./ gC2DInf{1}(:,r), 0, KInf, KInf ), PhiPhi2D(:,:,r));
  globOExp{2} = globOExp{2} + beta * kron( spdiags( sin(gInf.Beta).*(gInf.H+2*Q1(r)*gInf.M) .* gC2DInf{2}(:,r)  ./ gC2DInf{1}(:,r), 0, KInf, KInf ), PhiPhi2D(:,:,r));
end % for
end % function
