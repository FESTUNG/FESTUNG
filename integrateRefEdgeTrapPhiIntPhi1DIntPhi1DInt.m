function ret = integrateRefEdgeTrapPhiIntPhi1DIntPhi1DInt(N, qOrd, basesOnQuad2D, basesOnQuad1D)
[~, W] = quadRule1D(qOrd);
ret = zeros(N(1), N(2), N(2), 4);
for n = 1 : 2
  for i = 1 : N(1)
    for j = 1 : N(2)
      for l = 1 : N(2)
        ret(i,j,l,n) = W * ( basesOnQuad2D.phi1D(:,i,n) .* basesOnQuad1D.phi1D(:,j) .* basesOnQuad1D.phi1D(:,l) );
      end % for l
    end  % for j
  end  % for i
end  % for n
for n = 3 : 4
  for i = 1 : N(1)
    for j = 1 : N(2)
      for l = 1 : N(2)
        ret(i,j,l,n) = W * ( basesOnQuad2D.phi1D(:,i,n) .* basesOnQuad1D.phi0D(j,5-n) * basesOnQuad1D.phi0D(l,5-n) );
      end % for l
    end  % for j
  end  % for i
end  % for n
end  % function