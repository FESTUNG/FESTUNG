function ret = integrateRefEdge1DPhiIntPhiExtPhiExt(N, qOrd, basesOnQuad)
ret = { zeros(N, N, N, 4), zeros(N, N, N, 4) };
for n = 1 : 2
  for i = 1 : N
    for j = 1 : N
      for l = 1 : N
        ret{1}(i,j,l,n) = basesOnQuad.phi0D{qOrd}(i,n) * basesOnQuad.phi0D{qOrd}(j,3-n) * basesOnQuad.phi0D{qOrd}(l,3-n);
        ret{2}(i,j,l,n) = (2-n) * basesOnQuad.phi0D{qOrd}(i,n) * basesOnQuad.phi0D{qOrd}(j,3-n) * basesOnQuad.phi0D{qOrd}(l,3-n);
      end % for l
    end % for j
%     ret{1}(i,:,:,n) = (basesOnQuad.phi0D(i,n) * basesOnQuad.phi0D(:,3-n)) * basesOnQuad.phi0D(:,3-n).';
%     ret{2}(i,:,:,n) = ( ((n-1) * basesOnQuad.phi0D(i,n)) * basesOnQuad.phi0D(:,3-n) ) * basesOnQuad.phi0D(:,3-n).';
  end  % for i
end  % for n
end  % function