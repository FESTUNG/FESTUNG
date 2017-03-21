function ret = integrateRefEdgeTrapPhiIntPhiIntPerQuad(N, qOrd, basesOnQuad)
[~, W] = quadRule1D(qOrd); R = length(W);
ret = zeros(N, N, R, 4);
for n = 1 : 4
  for i = 1 : N
    for j = 1 : N
      ret(i,j,:,n) = W .* (basesOnQuad.phi1D{qOrd}(:,i,n) .* basesOnQuad.phi1D{qOrd}(:,j,n)).';
    end % for j
  end  % for i
end  % for n
end  % function