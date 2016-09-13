function globH = assembleMatElemTrapDphiPhi(g, refElemDphiPhi)
K = g.numT; N = size(refElemDphiPhi{1}, 1);
globH = { sparse(K*N, K*N), sparse(K*N, K*N) };
for m = 1 : 2
  for s = 1 : 3
    globH{m} = globH{m} + ...
               kron(spdiags(g.J0T{s}(:,3-m,3-m), 0, K, K), refElemDphiPhi{s}(:,:,  m)) - ...
               kron(spdiags(g.J0T{s}(:,3-m,  m), 0, K, K), refElemDphiPhi{s}(:,:,3-m));
  end % for s
end % for m
end  % function