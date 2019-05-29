function globR = assembleGlobRInf(gInf, hatRdiagInf, hatRoffdiagInf, cDGInf, gConst)
%% assembles global interior edge contributions of 0.5 * gConst * \int_E \varphi_i * H * \varphi_j * \nu
%% multiplying globR with H corresponds to interior edge integrals in weak form of 0.5 * gConst * H^2 in PDE
%% here the averages of the unknown height on the interior edges are used
[KInf, NInf, ~] = size(cDGInf);  globR = cell(2, 1);
globR{1} = sparse(KInf*NInf, KInf*NInf); globR{2} = sparse(KInf*NInf, KInf*NInf);
for nn = 1 : 3
  Rkn = 0.25 * gConst * gInf.detE0I(:,nn);
  for np = 1 : 3
    if nn == 2 || np == 2
      continue;
    end % if
    RknNu1 = bsxfun(@times, gInf.markE0IE0I{nn, np}, Rkn .* gInf.nuE0Inf(:,nn,1));
    RknNu2 = bsxfun(@times, gInf.markE0IE0I{nn, np}, Rkn .* gInf.nuE0Inf(:,nn,2));
    for l = 1 : NInf
      globR{1} = globR{1} + kron(bsxfun(@times, RknNu1, cDGInf(:,l,1).'), hatRoffdiagInf(:,:,l,nn,np));
      globR{2} = globR{2} + kron(bsxfun(@times, RknNu2, cDGInf(:,l,1).'), hatRoffdiagInf(:,:,l,nn,np));
    end % for
  end % for
  Rkn =  Rkn .* gInf.markE0Iint(:, nn);
  for l = 1 : NInf
    globR{1} = globR{1} + kron(spdiags(Rkn .* gInf.nuE0Inf(:,nn,1) .* cDGInf(:,l,1), 0, KInf, KInf), hatRdiagInf(:,:,l,nn));
    globR{2} = globR{2} + kron(spdiags(Rkn .* gInf.nuE0Inf(:,nn,2) .* cDGInf(:,l,1), 0, KInf, KInf), hatRdiagInf(:,:,l,nn));
  end % for
end % for
end % function
