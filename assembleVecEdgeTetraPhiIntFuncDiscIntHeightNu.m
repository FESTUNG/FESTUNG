function ret = assembleVecEdgeTetraPhiIntFuncDiscIntHeightNu(g2D, g1D, markE0T, dataDisc, heightV0T1D, N, qOrd, basesOnQuad)
K = g2D.numT;
[~, W] = quadRule1D(qOrd);
ret = zeros(K, N);
for n = 3 : 4
  funcQ0E = dataDisc * basesOnQuad.phi1D{qOrd}(:, :, n).';
  areaNuE0THeight = markE0T(:,n) .* g2D.areaE0T(:,n) .* g2D.nuE0T(:,n,1) ./ (g1D.markT2DT * heightV0T1D(:,5-n));
  ret = ret + bsxfun(@times, areaNuE0THeight, funcQ0E * (repmat(W(:), 1, N) .* basesOnQuad.phi1D{qOrd}(:, :, n)));
end  % for n
ret = reshape(ret.', K*N, 1);
end  % function