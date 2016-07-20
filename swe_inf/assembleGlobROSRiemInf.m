function globROSRiem = assembleGlobROSRiemInf(gInf, hatRdiag, cDGInf, gConst)
%% assembles global open sea boundary edge contributions of 0.5 * gConst * \int_E \varphi_i * H * \varphi_j * nu
%% multiplying globROSRiem with H corresponds to open sea boundary edge integrals in weak form of 0.5 * gConst * H^2 in PDE
%% here half of the values of the unknwon height from the interior are used (later they are used to form the average of interior and open sea boundary values)
[KInf, NInf, ~] = size(cDGInf);  globROSRiem = cell(2, 1);
globROSRiem{1} = sparse(KInf*NInf, KInf*NInf); globROSRiem{2} = sparse(KInf*NInf, KInf*NInf);
for n = 2 : 2
  Rkn = 0.25 * gConst * gInf.detE0I(:,n);
  for l = 1 : NInf
    globROSRiem{1} = globROSRiem{1} + kron(spdiags(Rkn .* gInf.nuE0Inf(:,n,1) .* cDGInf(:,l,1), 0, KInf, KInf), hatRdiag(:,:,l,n));
    globROSRiem{2} = globROSRiem{2} + kron(spdiags(Rkn .* gInf.nuE0Inf(:,n,2) .* cDGInf(:,l,1), 0, KInf, KInf), hatRdiag(:,:,l,n));
  end % for
end % for
end % function
