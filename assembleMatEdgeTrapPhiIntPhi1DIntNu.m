function ret = assembleMatEdgeTrapPhiIntPhi1DIntNu(g2D, g1D, markE0T, refEdgePhiIntPhi1DInt)
K = g2D.numT; barK = g1D.numT;
[N, barN, ~] = size(refEdgePhiIntPhi1DInt);
ret = { sparse(K*N, barK*barN), sparse(K*N, barK*barN) };
for n = 1 : 4
  areaE0Tint = markE0T(:,n) .* g2D.areaE0T(:,n);
  for m = 1 : 2
    areaNuE0Tint = areaE0Tint .* g2D.nuE0T(:,n,m);
    ret{m} = ret{m} + kron(bsxfun(@times, g1D.markT2DT, areaNuE0Tint), refEdgePhiIntPhi1DInt(:,:,n));
  end % for m
end  % for n
end  % function