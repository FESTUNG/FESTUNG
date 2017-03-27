function ret = assembleVecEdgeTrapPhiIntFuncCont1DNu(g, markE0T, funcCont, N, qOrd, basesOnQuad)
K = g.numT;
[Q, W] = quadRule1D(qOrd);
ret = { zeros(K, N), zeros(K, N) };
for n = 1 : 4
  [Q1, Q2] = gammaMapTetra(n, Q);
  funcQ0E = funcCont(g.mapRef2Phy(1, Q1, Q2));
  markAreaE0T = markE0T(:, n) .* g.areaE0T(:, n);
  for m = 1 : 2
    markAreaNuE0T = markAreaE0T .* g.nuE0T(:, n, m);
    for i = 1 : N
      ret{m}(:, i) = ret{m}(:, i) + markAreaNuE0T .* (funcQ0E * (W.' .* basesOnQuad.phi1D{qOrd}(:, i, n)) );
    end % for i
  end % for m
end % for n
ret = cellfun(@(c) reshape(c.', K*N, 1), ret, 'UniformOutput', false);
end % function