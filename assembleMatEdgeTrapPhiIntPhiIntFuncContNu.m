function ret = assembleMatEdgeTrapPhiIntPhiIntFuncContNu(g, markE0T, funcCont, qOrd, refEdgePhiIntPhiIntPerQuad)
K = g.numT;
[N, ~, ~, ~] = size(refEdgePhiIntPhiIntPerQuad);
[Q, ~] = quadRule1D(qOrd); R = length(Q);
ret = { sparse(K*N, K*N), sparse(K*N, K*N) };
for n = 1 : 4
  [Q1, Q2] = gammaMapQuad(n, Q);
  funcQ0E = funcCont(g.mapRef2Phy(1, Q1, Q2), g.mapRef2Phy(2, Q1, Q2));
  for m = 1 : 2
    markAreaNuE0T = markE0T(:, n) .* g.areaE0T(:, n) .* g.nuE0T(:, n, m);
    for r = 1 : R
      ret{m} = ret{m} + kron(spdiags(markAreaNuE0T .* funcQ0E(:, r), 0, K, K), refEdgePhiIntPhiIntPerQuad(:,:,r,n));
    end % for r
  end % for m
end % for n
end % function