function globQOS = assembleGlobQOSInf(gInf, hatQdiag)
%% assembles global open sea boundary edge contributions of \int_E \varphi_i \varphi_j * \nu
%% multiplying globQOS with uH, vH corresponds to opeb sea boundary edge integrals in weak form of div( [uH, vH]^T ) in PDE
%% here the values of the unknown momenta from the interior are used
K = gInf.numInfElem;  N = size(hatQdiag, 1);
markE0IOS = zeros(gInf.numInfElem,3);
markE0IOS(:,2) = 1;
globQOS = cell(2, 1);
globQOS{1} = sparse(K*N, K*N); globQOS{2} = sparse(K*N, K*N);
for n = 1 : 3
  QNkn = markE0IOS(:,n) .* gInf.detE0I(:,n);
  globQOS{1} = globQOS{1} + kron(spdiags(QNkn .* gInf.nuE0Inf(:,n,1), 0, K, K), hatQdiag(:,:,n));
  globQOS{2} = globQOS{2} + kron(spdiags(QNkn .* gInf.nuE0Inf(:,n,2), 0, K, K), hatQdiag(:,:,n));
end % for
end % function
