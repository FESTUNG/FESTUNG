function ret = integrateRefEdge1DPhiIntPhiIntPhiInt(N, qOrd, basesOnQuad)
ret = { zeros(N, N, N, 4), zeros(N, N, N, 4) };
for n = 1 : 2
  for i = 1 : N
    for j = 1 : N
      for l = 1 : N
      ret{1}(i,j,l,n) = basesOnQuad.phi0D{qOrd}(i,n) * basesOnQuad.phi0D{qOrd}(j,n) * basesOnQuad.phi0D{qOrd}(l,n);
      ret{2}(i,j,l,n) = (n-1) * basesOnQuad.phi0D{qOrd}(i,n) * basesOnQuad.phi0D{qOrd}(j,n) * basesOnQuad.phi0D{qOrd}(l,n);
      end % for l
    end % for j
%       ret{1}(i,:,:,n) = (basesOnQuad.phi0D(i,n) * basesOnQuad.phi0D(:,n)) * basesOnQuad.phi0D(:,n).';
%       ret{2}(i,:,:,n) = ( ((n-1) * basesOnQuad.phi0D(i,n)) * basesOnQuad.phi0D(:,n) ) * basesOnQuad.phi0D(:,n).';
  end  % for i
end  % for n
end  % function