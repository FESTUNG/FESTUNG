function hatG = integrateRefElemTrapDphiPhiPhi(N, qOrd, basesOnQuad)
[Q1, Q2, W] = gaussQuadRule2D(qOrd);
hatG = { zeros(N, N, N, 2), zeros(N, N, N, 2), zeros(N, N, N, 2) };
for i = 1 : N
  for j = 1 : N
    for l = 1 : j
      for m = 1 : 2
        ind = sub2ind([N, N, N, 2], [i i], [j l], [l j], [m m]);
        gradPhi2Dphi2Dphi2D = basesOnQuad.gradPhi2D(:,i,m) .* basesOnQuad.phi2D(:,j) .* basesOnQuad.phi2D(:,l);
        hatG{1}(ind) = W * gradPhi2Dphi2Dphi2D;
        hatG{2}(ind) = W * (gradPhi2Dphi2Dphi2D .* Q1);
        hatG{3}(ind) = W * (gradPhi2Dphi2Dphi2D .* Q2);
      end  % for m
    end  % for k
  end  % for j
end  % for i
end  % function