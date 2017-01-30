function ret = assembleMatEdgeTrapPhiIntPhi1DIntFuncDisc1DIntNuHeight(g2D, g1D, dataDisc1D, heightV0T1D, markE0Tint, refEdgePhiIntPhi1DIntPhi1DInt)
K = g2D.numT; barK = g1D.numT;
[N, ~, barN, ~] = size(refEdgePhiIntPhi1DIntPhi1DInt);
ret = sparse(barK*barN, K*N);
for n = 3 : 4
  areaNuE0Tint = markE0Tint(:,n) .* g2D.areaE0T(:,n) .* g2D.nuE0T(:,n,1);
  areaNuE0THeightint = areaNuE0Tint ./ (g1D.markT2DT * heightV0T1D(:,5-n));
  for l = 1 : barN
    areaNuE0THeightDataDisc = dataDisc1D(:,l) * areaNuE0THeightint.';
    ret = ret + kron(g1D.markT2DT.' .* areaNuE0THeightDataDisc, refEdgePhiIntPhi1DIntPhi1DInt(:,:,l,n).');
  end % for l
end  % for n
end  % function