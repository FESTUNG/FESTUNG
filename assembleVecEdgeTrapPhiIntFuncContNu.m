function ret = assembleVecEdgeTrapPhiIntFuncContNu(g, markE0Tbdr, funcCont, N, qOrd, basesOnQuad)
if iscell(funcCont)
  ret = assembleVecEdgeTrapPhiIntFuncContNuVector(g, markE0Tbdr, funcCont, N, qOrd, basesOnQuad);
else
  ret = assembleVecEdgeTrapPhiIntFuncContNuScalar(g, markE0Tbdr, funcCont, N, qOrd, basesOnQuad);
end % if
end % function

function ret = assembleVecEdgeTrapPhiIntFuncContNuVector(g, markE0Tbdr, funcCont, N, qOrd, basesOnQuad)
K = g.numT;
[Q, W] = quadRule1D(qOrd);
ret = { zeros(K, N), zeros(K, N) };
for n = 1 : 4
  [Q1, Q2] = gammaMapQuad(n, Q);
  funcQ0E = { funcCont{1}(g.mapRef2Phy(1, Q1, Q2), g.mapRef2Phy(2, Q1, Q2)), ...
              funcCont{2}(g.mapRef2Phy(1, Q1, Q2), g.mapRef2Phy(2, Q1, Q2)) };
  for m = 1 : 2
    markAreaNuE0T = markE0Tbdr(:, n) .* g.areaE0T(:, n) .* g.nuE0T(:, n, m);
    for i = 1 : N
      ret{m}(:, i) = ret{m}(:, i) + markAreaNuE0T .* ( funcQ0E{m} * (W.' .* basesOnQuad.phi1D{qOrd}(:, i, n)) );
    end % for i
  end  % for m
end  % for n
ret = cellfun(@(c) reshape(c.', K*N, 1), ret, 'UniformOutput', false);
end % function

function ret = assembleVecEdgeTrapPhiIntFuncContNuScalar(g, markE0Tbdr, funcCont, N, qOrd, basesOnQuad)
K = g.numT;
[Q, W] = quadRule1D(qOrd);
ret = { zeros(K, N), zeros(K, N) };
for n = 1 : 4
  [Q1, Q2] = gammaMapQuad(n, Q);
  funcQ0E = funcCont(g.mapRef2Phy(1, Q1, Q2), g.mapRef2Phy(2, Q1, Q2));
  for i = 1 : N
    markAreaIntE0T = markE0Tbdr(:, n) .* g.areaE0T(:, n) .* ( funcQ0E * (W.' .* basesOnQuad.phi1D{qOrd}(:, i, n)) );
    for m = 1 : 2
      ret{m}(:, i) = ret{m}(:, i) + markAreaIntE0T .* g.nuE0T(:, n, m);
    end % for m
  end  % for i
end  % for n
ret = cellfun(@(c) reshape(c.', K*N, 1), ret, 'UniformOutput', false);
end  % function