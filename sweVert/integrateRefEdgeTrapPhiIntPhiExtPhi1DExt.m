function ret = integrateRefEdgeTrapPhiIntPhiExtPhi1DExt(N, qOrd, basesOnQuad2D, basesOnQuad1D)
[~, W] = quadRule1D(qOrd);
mapE0E = [2 1 4 3];
ret = zeros(N(1), N(1), N(2), 4);
for n = 1 : 2
  np = mapE0E(n);
  for i = 1 : N(1)
    for j = 1 : N(1)
      for l = 1 : N(2)
        ret(i,j,l,n) = W * ( basesOnQuad2D.phi1D{qOrd}(:,i,n) .* basesOnQuad2D.phi1D{qOrd}(:,j,np) .* basesOnQuad1D.phi1D{qOrd}(:,l) );
      end % for l
    end  % for j
  end  % for i
end  % for n
for n = 3 : 4
  np = mapE0E(n);
  for i = 1 : N(1)
    for j = 1 : N(1)
      for l = 1 : N(2)
        ret(i,j,l,n) = W * ( basesOnQuad2D.phi1D{qOrd}(:,i,n) .* basesOnQuad2D.phi1D{qOrd}(:,j,np) * basesOnQuad1D.phi0D{qOrd}(l,n-2) );
      end % for l
    end  % for j
  end  % for i
end  % for n
end  % function