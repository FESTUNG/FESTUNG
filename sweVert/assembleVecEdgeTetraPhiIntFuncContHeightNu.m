function ret = assembleVecEdgeTetraPhiIntFuncContHeightNu(g2D, g1D, markE0T, funcCont, heightV0T1D, N, qOrd, basesOnQuad)
K = g2D.numT;
[Q, W] = quadRule1D(qOrd);
ret = zeros(K, N);
for n = 3 : 4
  [Q1, Q2] = gammaMapTetra(n, Q);
  funcQ0E = funcCont(g2D.mapRef2Phy(1, Q1, Q2), g2D.mapRef2Phy(2, Q1, Q2));
  areaNuE0T = markE0T(:,n) .* g2D.areaE0T(:,n) .* g2D.nuE0T(:,n,1);
  areaNuE0THeight = areaNuE0T ./ (g1D.markT2DT * heightV0T1D(:,5-n));
  for i = 1 : N
    ret(:, i) = ret(:, i) + areaNuE0THeight .* ( funcQ0E * (W.' .* basesOnQuad.phi1D{qOrd}(:, i, n)) );
  end  % for i
end  % for n
ret = reshape(ret.', K*N, 1);
end  % function