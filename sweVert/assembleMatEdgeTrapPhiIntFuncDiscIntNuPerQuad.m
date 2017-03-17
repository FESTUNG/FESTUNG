function ret = assembleMatEdgeTrapPhiIntFuncDiscIntNuPerQuad(g, markE0T, dataDisc, refEdgePhiIntPhiIntPerQuad)
[K, N] = size(dataDisc);
[~, ~, R, ~] = size(refEdgePhiIntPhiIntPerQuad);
ret = cell(4,2);
for n = 1 : 4
  markAreaE0T = markE0T(:, n) .* g.areaE0T(:, n);
  for m = 1 : 2
    ret{n, m} = sparse(K*N, K*R);
    markAreaNuE0T = markAreaE0T .* g.nuE0T(:, n, m);
    for j = 1 : N
      ret{n, m} = ret{n, m} + kron(spdiags(markAreaNuE0T .* dataDisc(:, j), 0, K, K), refEdgePhiIntPhiIntPerQuad(:,j,:,n));
    end % for j
  end % for m
end % for n
end % function