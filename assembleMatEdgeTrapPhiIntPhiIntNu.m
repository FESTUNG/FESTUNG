function ret = assembleMatEdgeTrapPhiIntPhiIntNu(g, markE0Tbdr, refEdgePhiIntPhiInt)
K = g.numT; N = size(refEdgePhiIntPhiInt, 1);
ret = { sparse(K*N, K*N), sparse(K*N, K*N) };
for n = 1 : 4
  areaE0Tbdr = markE0Tbdr(:,n) .* g.areaE0T(:,n);
  for m = 1 : 2
    ret{m} = ret{m} + kron(spdiags(areaE0Tbdr .* g.nuE0T(:,n,m), 0, K, K), refEdgePhiIntPhiInt(:,:,n));
  end % for m
end % for n
end  % function