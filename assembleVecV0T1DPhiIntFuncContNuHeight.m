function ret = assembleVecV0T1DPhiIntFuncContNuHeight(g2D, g1D, markE0T, funcCont, heightV0T1D, barN, qOrd, basesOnQuad)
barK = g1D.numT;
[Q, W] = quadRule1D(qOrd);
ret = zeros(barK, barN);
for n = 3 : 4
  [Q1, Q2] = gammaMapTetra(n, Q);
  funcQ0E = funcCont(g2D.mapRef2Phy(1, Q1, Q2), g2D.mapRef2Phy(2, Q1, Q2));
  areaNuE0T = markE0T(:,n) .* g2D.areaE0T(:,n) .* g2D.nuE0T(:,n,1);
  integrand = g1D.markT2DT.' * (areaNuE0T .* (funcQ0E * W.'));
  ret = ret + (integrand ./ heightV0T1D(:, 5-n)) * basesOnQuad.phi0D{qOrd}(:, 5-n).';
end  % for n
ret = reshape(ret.', barK*barN, 1);
end  % function