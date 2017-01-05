function ret = assembleMatElem1DDphiPhiFuncDisc(g, dataDisc, refElemDphiPhiPhi)
[K, N] = size(dataDisc{1});
ret = sparse(K*N, K*N);
for s = 1 : 2
  for l = 1 : N
    ret = ret + kron(spdiags(dataDisc{s}(:,l) .* g.detJ0T, 0, K, K), refElemDphiPhiPhi{s}(:,:,l));
  end % for l
end % for s
end % function