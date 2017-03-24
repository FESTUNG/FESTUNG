function ret = integrateRefEdgeTrapPhiIntPhi1DExt(N, qOrd, basesOnQuad2D, basesOnQuad1D)
[~, W] = quadRule1D(qOrd);
ret = zeros(N(1), N(2), 4);
for n = 1 : 2
  for i = 1 : N(1)
    for j = 1 : N(2)
      ret(i,j,n) = W * ( basesOnQuad2D.phi1D{qOrd}(:,i,n) .* basesOnQuad1D.phi1D{qOrd}(:,j) );
    end  % for j
  end  % for i
end  % for n
for n = 3 : 4
  for i = 1 : N(1)
    for j = 1 : N(2)
      ret(i,j,n) = W * ( basesOnQuad2D.phi1D{qOrd}(:,i,n) * basesOnQuad1D.phi0D{qOrd}(j,n-2) );
    end  % for j
  end  % for i
end  % for n
end  % function