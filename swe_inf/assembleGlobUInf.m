function globU = assembleGlobUInf(gInf, PhiPhi1D, PhiPhioffdiag, cDG, p, beta, gCdiagInf, gCoffdiagInf)
%% assembles global interior edge contributions of \int_E \varphi_i * ( uH * nu^{1} + vH * nu^{2} ) / H * \varphi_j
%% multiplying globU with uH (in the first momentum equation), vH (in the second momentum equation)
%% corresponds to interior edge integrals in weak form of div( [uH*uH, vH*uH; uH*vH, vH*vH] / H ) in PDE
%% here the averages of the unknown momenta on the interior edges are used
% global gCdiagInf gCoffdiagInf
[KInf, NInf, ~] = size(cDG);
globU = sparse(KInf*NInf, KInf*NInf);
pInf = NInf/(p+1)-1;  qOrd = max(2*pInf+1, 1); [~, W] = quadRuleInf1D(qOrd, beta);
for nn = 1 : 3
  Ukn = 0.5 * gInf.detE0I(:,nn);
  for np = 1 : 3
    if nn == 2 || np == 2
      continue;
    end
    u = gCoffdiagInf{2,nn,np}(:,:) ./ gCoffdiagInf{1,nn,np}(:,:);
    u = max(u,0) + min(u,0);
    v = gCoffdiagInf{3,nn,np}(:,:) ./ gCoffdiagInf{1,nn,np}(:,:);
    v = max(v,0) + min(v,0);
    for r = 1 : length(W)
      globU = globU + kron( bsxfun( @times, gInf.markE0IE0I{nn, np}, Ukn .* gInf.markE0Iaux{nn,np} .* ( gInf.nuE0Inf(:,nn,1) .* u(:,r) + gInf.nuE0Inf(:,nn,2) .* v(:,r) ) ), PhiPhioffdiag(:,:,r,nn,np) );
    end % for
  end % for
  for r = 1 : length(W)
    globU = globU + kron( spdiags( Ukn .* gInf.markE0Iint(:,nn) .* ( gInf.nuE0Inf(:,nn,1) .* gCdiagInf{2,nn}(:,r) + gInf.nuE0Inf(:,nn,2) .* gCdiagInf{3,nn}(:,r) ) ./ gCdiagInf{1,nn}(:,r), 0, KInf, KInf), PhiPhi1D(:,:,r,nn) );
  end % for
end % for
end % function
