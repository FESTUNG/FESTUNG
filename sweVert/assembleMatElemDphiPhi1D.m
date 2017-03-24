function ret = assembleMatElemDphiPhi1D(g, refElemDphiPhi1D)
K = g.numT; barK = g.g1D.numT; N = size(refElemDphiPhi1D{1});
ret = { sparse(K*N(1), barK*N(2)), sparse(K*N(1), barK*N(2)) };
for m = 1 : 2
  for s = 1 : 3
    ret{m} = ret{m} + ...
             kron(bsxfun(@times, g.g1D.markT2DT, g.J0T{s}(:,3-m,3-m)), refElemDphiPhi1D{s}(:,:,  m)) - ...
             kron(bsxfun(@times, g.g1D.markT2DT, g.J0T{s}(:,3-m,  m)), refElemDphiPhi1D{s}(:,:,3-m));
  end % for s
end % for m
end % function