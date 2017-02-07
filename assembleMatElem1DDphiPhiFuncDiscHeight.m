function ret = assembleMatElem1DDphiPhiFuncDiscHeight(dataDisc, heightQ0T, refElemDphiPhiPhiPerQuad)
[K, N] = size(dataDisc{1}); R = size(heightQ0T, 2);
ret = sparse(K*N, K*N);
invHeightQ0T = 1./heightQ0T;
for s = 1 : 2
  for l = 1 : N
    refElemDphiPhiPhiHeight = zeros(K*N, N);
    for r = 1 : R
      refElemDphiPhiPhiHeight = refElemDphiPhiPhiHeight + kron(invHeightQ0T(:,r), refElemDphiPhiPhiPerQuad{s}(:,:,l,r));
    end % for r
    ret = ret + kronVec(spdiags(dataDisc{s}(:,l), 0, K, K), refElemDphiPhiPhiHeight);
  end % for l
end % for s
end % function