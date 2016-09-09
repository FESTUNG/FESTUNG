function hatH = integrateRefElemTrapDphiPhi(N, qOrd, basesOnQuad)
[Q1, Q2, W] = gaussQuadRule2D(qOrd);
hatH = { zeros(N, N, 2), zeros(N, N, 2), zeros(N, N, 2) };
for i = 1 : N
  for j = 1 : N
    for m = 1 : 2
      phi2DgradPhi2D = basesOnQuad.phi2D(:,j) .* basesOnQuad.gradPhi2D(:,i,m);
      hatH{1}(i,j,m) = W * phi2DgradPhi2D;
      hatH{2}(i,j,m) = W * ( phi2DgradPhi2D .* Q1 );
      hatH{3}(i,j,m) = W * ( phi2DgradPhi2D .* Q2 );
    end  % for m
  end  % for j
end  % for i
end  % function