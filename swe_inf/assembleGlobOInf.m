function globO = assembleGlobOInf(gInf, gradPhiPhi2D, cDG, p, beta, gC2DInf)
%% assembles global contributions of - \int_T ( ( \partial_{x^1} \varphi_i ) * uH / H + ( \partial_{x^2} \varphi_i ) * vH / H ) * \varphi_j
%% multiplying globO with uH (in the first momentum equation), vH (in the second momentum equation)
%% corresponds to area integrals in weak form of div( [uH*uH, uH*vH; vH*uH, vH*vH] / H ) in PDE
% global gC2DInf
[KInf, NInf, ~] = size(cDG);
qOrd = max(2*p,1);
pInf = NInf/(p+1)-1;  qOrdInf = max(2*pInf, 1); [Q1, Q2, W] = quadRuleInf2D(qOrd, qOrdInf, beta);
globO = sparse(KInf*NInf, KInf*NInf);
for r = 1 : length(W)
  globO = globO - ( kron( spdiags( ( cos(gInf.Beta).*(gInf.H+2*Q1(r)*gInf.M) .* gC2DInf{2}(:,r) + sin(gInf.Beta).*(gInf.H+2*Q1(r)*gInf.M) .* gC2DInf{3}(:,r) ) ./ gC2DInf{1}(:,r), 0, KInf, KInf ), gradPhiPhi2D(:,:,r,1) ) ...
                  + kron( spdiags( ( (-sin(gInf.Alpha)-cos(gInf.Beta)*2.*gInf.M*(Q2(r)-0.5)) .* gC2DInf{2}(:,r) + (cos(gInf.Alpha)-sin(gInf.Beta)*2.*gInf.M*(Q2(r)-0.5)) .* gC2DInf{3}(:,r) ) ./ gC2DInf{1}(:,r), 0, KInf, KInf ), gradPhiPhi2D(:,:,r,2) ) );
end % for
end % function
