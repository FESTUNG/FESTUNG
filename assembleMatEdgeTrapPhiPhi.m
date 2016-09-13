function ret = assembleMatEdgeTrapPhiPhi(g, markE0Tint, refEdgePhiIntPhiInt, refEdgePhiIntPhiExt)
K = g.numT; N = size(refEdgePhiIntPhiInt, 1);
ret = sparse(K*N, K*N);
for n = 1 : 4
  ret = ret + ...
        kron(spdiags(markE0Tint(:,n), 0, K, K), refEdgePhiIntPhiInt(:,:,n)) - ...
        kron(g.markE0TE0T{n}, refEdgePhiIntPhiExt(:,:,n));
end % for n
end % function