function globQ = assembleGlobQInf(gInf, hatQdiag, hatQoffdiag)
%% assembles global interior edge contributions of \int_E \varphi_i * \varphi_j * \nu
%% multiplying globQ with uH, vH corresponds to interior edge integrals in weak form of div( [uH, vH]^T ) in PDE
%% here the averages of the unknown momenta on the interior edges are used
NInf = size(hatQdiag,1); KInf = gInf.numInfElem;
globQ = cell(2, 1);
globQ{1} = sparse(KInf*NInf, KInf*NInf);  globQ{2} = sparse(KInf*NInf, KInf*NInf);
for nn = 1 : 3
  Qkn = 0.5 * gInf.detE0I(:,nn);
  for np = 1 : 3
    globQ{1} = globQ{1} + kron(bsxfun(@times, gInf.markE0IE0I{nn,np}, Qkn .* gInf.nuE0Inf(:,nn,1)), hatQoffdiag(:,:,nn,np));
    globQ{2} = globQ{2} + kron(bsxfun(@times, gInf.markE0IE0I{nn,np}, Qkn .* gInf.nuE0Inf(:,nn,2)), hatQoffdiag(:,:,nn,np));
  end % for
  Qkn = Qkn .* gInf.markE0Iint(:,nn);
  globQ{1} = globQ{1} + kron(spdiags(Qkn .* gInf.nuE0Inf(:,nn,1), 0, KInf, KInf), hatQdiag(:,:,nn));
  globQ{2} = globQ{2} + kron(spdiags(Qkn .* gInf.nuE0Inf(:,nn,2), 0, KInf, KInf), hatQdiag(:,:,nn));
end % for
end % function
