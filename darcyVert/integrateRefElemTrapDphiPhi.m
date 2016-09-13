function hatH = integrateRefElemTrapDphiPhi(N, qOrd, basesOnQuad)
[Q, W] = quadRule1D(qOrd); [Q1, Q2] = meshgrid(Q); W = W' * W;
Q1 = Q1(:); Q2 = Q2(:); W = W(:)';
hatH = { zeros(N, N, 2), zeros(N, N, 2), zeros(N, N, 2) };
for m = 1 : 2
  for j = 1 : N
    for i = 1 : N
      phi2DgradPhi2D = basesOnQuad.phi2D(:,j) .* basesOnQuad.gradPhi2D(:,i,m);
      hatH{1}(i,j,m) = W * phi2DgradPhi2D;
      hatH{2}(i,j,m) = W * ( phi2DgradPhi2D .* Q1 );
      hatH{3}(i,j,m) = W * ( phi2DgradPhi2D .* Q2 );
    end  % for i
  end  % for j
end  % for m
end  % function