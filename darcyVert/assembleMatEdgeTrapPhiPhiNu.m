function ret = assembleMatEdgeTrapPhiPhiNu(g, markE0Tint, refEdgePhiIntPhiInt, refEdgePhiIntPhiExt)
K = g.numT; N = size(refEdgePhiIntPhiInt, 1);
ret = { sparse(K*N, K*N), sparse(K*N, K*N) };
for n = 1 : 4
  areaE0Tint = 0.5 * markE0Tint(:,n) .* g.areaE0T(:,n);
  for m = 1 : 2
    areaNuE0Tint = areaE0Tint .* g.nuE0T(:,n,m);
    ret{m} = ret{m} + ...
             kron(spdiags(areaNuE0Tint, 0, K, K ) , refEdgePhiIntPhiInt(:,:,n)) + ...
             kron(bsxfun(@times, g.markE0TE0T{n}, areaNuE0Tint), refEdgePhiIntPhiExt(:,:,n));
  end % for m
end  % for n
end  % function