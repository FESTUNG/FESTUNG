function ret = assembleMatEdgeTrapPhiIntPhiInt(g, markE0Tbdr, refEdgePhiIntPhiInt, coefE0T)
if nargin < 4
  coefE0T = g.areaE0T;
end % if
K = g.numT; N = size(refEdgePhiIntPhiInt, 1);
ret = sparse(K*N, K*N);
for n = 1 : 4
  ret = ret + kron(spdiags(markE0Tbdr(:,n) .* coefE0T(:,n), 0, K, K), refEdgePhiIntPhiInt(:,:,n));
end  % for
end  % function