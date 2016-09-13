function ret = integrateRefEdgeTrapPhiIntPhiExtPhiExt(N, qOrd, basesOnQuad)
[~, W] = quadRule1D(qOrd);
mapE0E = [2 1 4 3];
ret = zeros(N, N, N, 4);
for n = 1 : 4
  for l = 1 : N
    for j = 1 : l
      for i = 1 : N
        ind = sub2ind([N N N 4], [i i], [j l], [l j], [n n]);
        ret(ind) = W * ( basesOnQuad.phi1D(:,i,n) .* basesOnQuad.phi1D(:,l,mapE0E(n)) .* basesOnQuad.phi1D(:,j,mapE0E(n)) );
      end  % for i
    end  % for j
  end  % for l
end % for n
end % function