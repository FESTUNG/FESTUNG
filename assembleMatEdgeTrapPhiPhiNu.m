function ret = assembleMatEdgeTrapPhiPhiNu(g, markE0Tint, refEdgePhiIntPhiInt, refEdgePhiIntPhiExt, edgeIndex, normalIndex)
if nargin < 5
  edgeIndex = 1 : 4;
end % if
if nargin < 6
  normalIndex = 1 : 2;
end % if
K = g.numT; N = size(refEdgePhiIntPhiInt, 1);
ret = cellfun(@(c) sparse(K*N, K*N), cell(size(normalIndex)), 'UniformOutput', false); %{ sparse(K*N, K*N), sparse(K*N, K*N) };
for n = edgeIndex
  areaE0Tint = 0.5 * markE0Tint(:,n) .* g.areaE0T(:,n);
  for m = normalIndex
    areaNuE0Tint = areaE0Tint .* g.nuE0T(:,n,m);
    ret{m} = ret{m} + ...
             kron(spdiags(areaNuE0Tint, 0, K, K ) , refEdgePhiIntPhiInt(:,:,n)) + ...
             kron(bsxfun(@times, g.markE0TE0T{n}, areaNuE0Tint), refEdgePhiIntPhiExt(:,:,n));
  end % for m
end  % for n
end  % function