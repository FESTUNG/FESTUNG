function ret = integrateRefEdgeTrapPhiIntPhiExt(N, qOrd, basesOnQuad)
[~, W] = quadRule1D(qOrd);
mapE0E = [2 1 4 3];
ret = zeros(N, N, 4);
for n = 1 : 4
  np = mapE0E(n);
  for j = 1 : N
    for i = 1 : N
      ret(i,j,n) = W * (basesOnQuad.phi1D{qOrd}(:,i,n) .* basesOnQuad.phi1D{qOrd}(:,j,np));
    end  % for j
  end  % for i
end % for nn
end  % function