function ret = assembleVecEdgeTrapPhiIntFuncDiscIntFuncContNu(g, markE0Tbdr, dataDisc, funcCont, refEdgePhiIntPhiIntPerQuad)
[K, N] = size(dataDisc);
[Q, W] = quadRule1D(qOrd);
ret = { zeros(K, N), zeros(K, N) };
for n = 1 : 4
  [Q1, Q2] = gammaMapQuad(n, Q);
  funcQ0E = funcCont(g.mapRef2Phy(1, Q1, Q2), g.mapRef2Phy(2, Q1, Q2));
  for i = 1 : N
    markAreaIntE0T = markE0Tbdr(:, n) .* g.areaE0T(:, n) .* ( funcQ0E * (W.' .* basesOnQuad.phi1D(:, i, n)) );
    for m = 1 : 2
      markAreaNuIntE0T = markAreaIntE0T .* g.nuE0T(:, n, m);
      for l = 1 : N
        ret{m}(:, i) = ret{m}(:, i) + markAreaIntE0T .* g.nuE0T(:, n, m);
      end % for l
    end % for m
  end  % for i
end  % for n
ret = cellfun(@(c) reshape(c.', K*N, 1), ret, 'UniformOutput', false);
end  % function