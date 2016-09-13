function ret = assembleMatEdgeTrapPhiPhiFuncDiscNu(g, markE0Tint, refEdgePhiIntPhiIntPhiInt, refEdgePhiIntPhiExtPhiExt, dataDisc)
[K, N] = size(dataDisc{1,1});
ret = { sparse(K*N, K*N), sparse(K*N, K*N) };
for n = 1 : 4
  for r = 1 : 2
    areaNuE0Tint = 0.5 * g.areaE0T(:,n) .* g.nuE0T(:,n,r) .* markE0Tint(:,n);
    for m = 1 : 2
      for l = 1 : N
        ret{m} = ret{m} + kron(spdiags(areaNuE0Tint .* dataDisc{r,m}(:,l), 0, K, K), ...
                               refEdgePhiIntPhiIntPhiInt(:,:,l,n)) + ...
                          kron(g.markE0TE0T{n} .* (areaNuE0Tint * dataDisc{r,m}(:,l)'), ...
                               refEdgePhiIntPhiExtPhiExt(:,:,l,n));
      end % for l
    end  % for m
  end % for r
end  % for n
end  % function