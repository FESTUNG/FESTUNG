function globROSRiem = assembleGlobROSRiem(g, markE0TOS, hatRdiag, cDG, gConst)
%% assembles global open sea boundary edge contributions of 0.5 * gConst * \int_E \varphi_i * H * \varphi_j * nu
%% multiplying globROSRiem with H corresponds to open sea boundary edge integrals in weak form of 0.5 * gConst * H^2 in PDE
%% here half of the values of the unknwon height from the interior are used (later they are used to form the average of interior and open sea boundary values)
[K, N, ~] = size(cDG);  globROSRiem = cell(2, 1);
globROSRiem{1} = sparse(K*N, K*N); globROSRiem{2} = sparse(K*N, K*N);
for n = 1 : 3
  Rkn = 0.25 * gConst * markE0TOS(:,n) .* g.areaE0T(:,n);
  for l = 1 : N
    globROSRiem{1} = globROSRiem{1} + kron(spdiags(Rkn .* g.nuE0T(:,n,1) .* cDG(:,l,1), 0, K, K), hatRdiag(:,:,l,n));
    globROSRiem{2} = globROSRiem{2} + kron(spdiags(Rkn .* g.nuE0T(:,n,2) .* cDG(:,l,1), 0, K, K), hatRdiag(:,:,l,n));
  end % for
end % for
end % function
