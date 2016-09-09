function hatM = integrateRefElemTrapPhiPhi(N, qOrd, basesOnQuad)
[Q1, Q2, W] = gaussQuadRule2D(qOrd);
hatM = { zeros(N, N), zeros(N, N), zeros(N, N) };
for i = 1 : N
  for j = 1 : i
    ind = sub2ind([N N], [i j], [j i]);
    phi2Dphi2D = basesOnQuad.phi2D(:, i) .* basesOnQuad.phi2D(:, j);
    hatM{1}(ind) = W * phi2Dphi2D;
    hatM{2} = W * ( phi2Dphi2D .* Q1 );
    hatM{3} = W * ( phi2Dphi2D .* Q2 );
  end  % for j
end  % for i
end  % function