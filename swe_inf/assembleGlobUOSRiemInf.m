function globUOSRiem = assembleGlobUOSRiemInf(gInf, PhiPhi1D, cDGInf, p, gCdiagInf)
%% assembles global open sea boundary edge contributions of \int_E \varphi_i * ( uH * nu^{1} + vH * nu^{2} ) / H * \varphi_j
%% multiplying globUOSRiem with uH (in the first momentum equation), vH (in the second momentum equation)
%% corresponds to open sea boundary edge integrals in weak form of div( [uH*uH, vH*uH; uH*vH, vH*vH] / H ) in PDE
%% here the values of the unknown height and momenta from the interior are used
% global gCdiagInf
[KInf, NInf, ~] = size(cDGInf);  globUOSRiem = sparse(KInf*NInf, KInf*NInf); 
qOrd = max(2*p+1, 1); [~, W] = quadRule1D(qOrd);
for n = 2 : 2
  UOSkn = 0.5 * gInf.detE0I(:,n);
  for r = 1 : length(W)
    globUOSRiem = globUOSRiem + kron( spdiags( UOSkn .* ( gCdiagInf{2,n}(:,r) .* gInf.nuE0Inf(:,n,1) + gCdiagInf{3,n}(:,r) .* gInf.nuE0Inf(:,n,2) ) ./ gCdiagInf{1,n}(:,r), 0, KInf, KInf ), PhiPhi1D(:,:,r,n) );
  end % for
end % for
end % function
