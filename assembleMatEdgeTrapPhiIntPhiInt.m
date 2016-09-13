function ret = assembleMatEdgeTrapPhiIntPhiInt(g, markE0Tbdr, refEdgePhiIntPhiInt)
K = g.numT; N = size(refEdgePhiIntPhiInt, 1);
ret = sparse(K*N, K*N);
for n = 1 : 4
  ret = ret + kron(spdiags(markE0Tbdr(:,n), 0, K, K), refEdgePhiIntPhiInt(:,:,n));
end  % for
end  % function