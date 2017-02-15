function ret = integrateRefEdgeTrapPhiIntPerQuad(N, qOrd, basesOnQuad)
[~, W] = quadRule1D(qOrd); R = length(W);
ret = zeros(N, R, 4);
for n = 1 : 4
  for i = 1 : N
    ret(i,:,n) = W .* basesOnQuad.phi1D(:,i,n).';
  end  % for i
end  % for n
end  % function