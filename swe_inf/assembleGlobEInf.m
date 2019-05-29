function globE = assembleGlobEInf(gInf, PhiPhi2D, cDG, tauConst, p, beta, gC2DInf)
%% assembles global contributions of tauConst * \int_T \varphi_i * (uH^2 + vH^2)^0.5 / H * \varphi_j
%% multiplying globE with uH, vH corresponds to area integrals in weak form of \tau_{bf} * uH, \tau_{bf} * vH respectively in PDE
% global gC2DInf
[KInf, NInf, ~] = size(cDG);
qOrd = max(2*p,1);
pInf = NInf/(p+1)-1; qOrdInf = max(2*pInf, 1); [~, ~, W] = quadRuleInf2D(qOrd, qOrdInf, beta);
globE = sparse(KInf*NInf, KInf*NInf);
for r = 1:length(W)
  globE = globE + tauConst * kron( spdiags( gInf.detQuadPoints(:,r) .* (gC2DInf{2}(:,r).^2  + gC2DInf{3}(:,r).^2).^0.5 ./ (gC2DInf{1}(:,r).^2) , 0, KInf, KInf ), PhiPhi2D(:,:,r) );
end % for
end % function
