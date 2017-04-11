function ret = assembleMatV0T1DPhiPhiFuncDiscNuHeight(g, dataDisc, heightV0T, markV0T, refEdgePhiIntPhiIntPhiInt, refEdgePhiIntPhiExtPhiExt)
[K,N] = size(dataDisc{1});
ret = sparse(K*N, K*N);
for n = 1 : 2
  nuV0THeight = 0.5 * markV0T(:,n) .* g.nuV0T(:,n) ./ heightV0T(:,n);
  for s = 1 : 2
    for l = 1 : N
      ret = ret + kron(spdiags(nuV0THeight .* dataDisc{s}(:,l), 0, K, K), ...
                       refEdgePhiIntPhiIntPhiInt{s}(:,:,l,n)) + ...
                  kron(g.markV0TV0T{n} .* (nuV0THeight * dataDisc{s}(:,l).'), ...
                       refEdgePhiIntPhiExtPhiExt{s}(:,:,l,n));
    end % for l
  end % for s
end % for n
end % function