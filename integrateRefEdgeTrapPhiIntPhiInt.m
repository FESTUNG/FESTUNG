function hatSdiag = integrateRefEdgeTrapPhiIntPhiInt(N, qOrd, basesOnQuad)
[~, W] = quadRule1D(qOrd);
hatSdiag = zeros(N, N, 4);
for n = 1 : 4
  for i = 1 : N
    for j = 1 : i
      ind = sub2ind([N N 4], [i j], [j i], [n n]);
      hatSdiag(ind) = W * ( basesOnQuad.phi1D(:,i,n) .* basesOnQuad.phi1D(:,j,n) );
    end  % for j
  end  % for i
end  % for n
end  % function