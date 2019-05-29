function globRL = assembleGlobRLInf(gInf, hatRdiag, cDGInf, gConst)
%% assembles global land boundary edge contributions of 0.5 * gConst * \int_E \varphi_i * H * \varphi_j * \nu
%% multiplying globRL with H corresponds to land boundary edge integrals in weak form of 0.5 * gConst * H^2 in PDE
%% here the values of the unknown height from the interior are used
[KInf, NInf, ~] = size(cDGInf);  globRL = cell(2, 1);
globRL{1} = sparse(KInf*NInf, KInf*NInf); globRL{2} = sparse(KInf*NInf, KInf*NInf);
for n = 1 : 3
  RLkn = gInf.markE0IL(:,n) .* gInf.detE0I(:,n);
  for l = 1 : NInf
    globRL{1} = globRL{1} + 0.5 * gConst * kron(spdiags(RLkn .* gInf.nuE0Inf(:,n,1) .* cDGInf(:,l,1), 0, KInf, KInf), hatRdiag(:,:,l,n));
    globRL{2} = globRL{2} + 0.5 * gConst * kron(spdiags(RLkn .* gInf.nuE0Inf(:,n,2) .* cDGInf(:,l,1), 0, KInf, KInf), hatRdiag(:,:,l,n));
  end % for
end % for
end % function
