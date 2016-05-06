function ret = assembleMatEdgePhiIntPerQuad(g, markE0Tbdr, refEdgePhiIntPerQuad, areaE0Tbdr)
K = g.numT;
ret = cell(3,3);
for n = 1 : 3
  if nargin > 3
    ret{n} = 0.5 * kron( spdiags(areaE0Tbdr{n}, 0, K, K), refEdgePhiIntPerQuad(:,:,n) );
  else
    ret{n} = 0.5 * kron( spdiags(markE0Tbdr(:,n) .* g.areaE0T(:,n), 0, K, K), refEdgePhiIntPerQuad(:,:,n) );
  end % if
end % for
end % function
