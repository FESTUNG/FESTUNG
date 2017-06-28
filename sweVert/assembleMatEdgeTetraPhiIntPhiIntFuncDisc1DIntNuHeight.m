function ret = assembleMatEdgeTetraPhiIntPhiIntFuncDisc1DIntNuHeight(g2D, g1D, dataDisc1D, heightV0T1D, markE0T, refEdgePhiIntPhiIntPhi1DInt)
K = g2D.numT;
[N, ~, barN, ~] = size(refEdgePhiIntPhiIntPhi1DInt);
ret = sparse(K*N, K*N);
dataDisc2D = g1D.markT2DT * dataDisc1D;
for n = 3 : 4
  areaNuE0Tint = markE0T(:,n) .* g2D.areaE0T(:,n) .* g2D.nuE0T(:,n,1);
  areaNuE0THeightint = areaNuE0Tint ./ (g1D.markT2DT * heightV0T1D(:,5-n));
  for l = 1 : barN
    ret = ret + kron(spdiags(areaNuE0THeightint .* dataDisc2D(:,l), 0, K, K), refEdgePhiIntPhiIntPhi1DInt(:,:,l,n));
  end % for l
end  % for n
end  % function