function ret = assembleMatEdgeTrapPhiPhi(g, markE0Tint, refEdgePhiIntPhiInt, refEdgePhiIntPhiExt, coefE0T)
if nargin < 5
  coefE0T = g.areaE0T;
end % if
K = g.numT; N = size(refEdgePhiIntPhiInt, 1);
ret = sparse(K*N, K*N);
for n = 1 : 4
  areaE0Tint = markE0Tint(:,n) .* coefE0T(:,n);
  ret = ret + ...
        kron(spdiags(areaE0Tint, 0, K, K), refEdgePhiIntPhiInt(:,:,n)) - ...
        kron(bsxfun(@times, g.markE0TE0T{n}, areaE0Tint), refEdgePhiIntPhiExt(:,:,n));
end % for n
end % function