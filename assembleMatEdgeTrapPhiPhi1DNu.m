function ret = assembleMatEdgeTrapPhiPhi1DNu(g2D, g1D, markE0Tint, refEdgePhiIntPhi1DInt, refEdgePhiIntPhi1DExt)
K = g2D.numT; barK = g1D.numT;
[N, barN, ~] = size(refEdgePhiIntPhi1DInt);
ret = { sparse(K*N, barK*barN), sparse(K*N, barK*barN) };
for n = 1 : 4
  areaE0Tint = 0.5 * markE0Tint(:,n) .* g2D.areaE0T(:,n);
  for m = 1 : 2
    areaNuE0Tint = areaE0Tint .* g2D.nuE0T(:,n,m);
    ret{m} = ret{m} + ...
             kron(bsxfun(@times, g1D.markT2DT, areaNuE0Tint), refEdgePhiIntPhi1DInt(:,:,n)) + ...
             kron(bsxfun(@times, g2D.markE0TE0T{n} * double(g1D.markT2DT), areaNuE0Tint), refEdgePhiIntPhi1DExt(:,:,n));
  end % for m
end  % for n
end  % function