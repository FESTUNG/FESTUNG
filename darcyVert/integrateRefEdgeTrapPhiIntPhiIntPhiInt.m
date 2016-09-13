function ret = integrateRefEdgeTrapPhiIntPhiIntPhiInt(N, qOrd, basesOnQuad)
[~, W] = quadRule1D(qOrd);
ret = zeros(N, N, N, 4);
for n = 1 : 4  
  for l = 1 : N
    for i = 1 : l
      for j = 1 : i
        ind = sub2ind([N N N 4], [i j i j l l], [j i l l i j], [l l j i j i], [n n n n n n]);
        ret(ind) = W * ( basesOnQuad.phi1D(:,i,n) .* basesOnQuad.phi1D(:,l,n) .* basesOnQuad.phi1D(:,j,n) );
      end  % for j
    end  % for i
  end  % for l
end  % for n
end  % function