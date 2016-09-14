function ret = assembleVecEdgeTrapPhiIntFuncCont(g, markAreaE0Tbdr, funcCont, N, qOrd, basesOnQuad)
K = g.numT;
[Q, W] = quadRule1D(qOrd);
ret = zeros(K, N);
for n = 1 : 4
  [Q1, Q2] = gammaMapTrap(n, Q);
  funcQ0E = funcCont(g.mapRef2Phy(1, Q1, Q2), g.mapRef2Phy(2, Q1, Q2));
  for i = 1 : N
    ret(:, i) = ret(:, i) + markAreaE0Tbdr(:, n) .* ( funcQ0E * ( W.' .* basesOnQuad.phi1D(:, i, n) ) );
  end  % for i
end  % for n
ret = reshape(ret.', K*N, 1);
end  % function