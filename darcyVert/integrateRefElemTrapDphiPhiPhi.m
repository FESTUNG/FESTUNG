function hatG = integrateRefElemTrapDphiPhiPhi(N, qOrd, basesOnQuad)
[Q, W] = quadRule1D(qOrd); [Q1, Q2] = meshgrid(Q); W = W' * W;
Q1 = Q1(:); Q2 = Q2(:); W = W(:)';
hatG = { zeros(N, N, N, 2), zeros(N, N, N, 2), zeros(N, N, N, 2) };
for m = 1 : 2
  for j = 1 : N
    for l = 1 : j
      for i = 1 : N
        ind = sub2ind([N, N, N, 2], [i i], [j l], [l j], [m m]);
        gradPhi2Dphi2Dphi2D = basesOnQuad.gradPhi2D(:,i,m) .* basesOnQuad.phi2D(:,j) .* basesOnQuad.phi2D(:,l);
        hatG{1}(ind) = W * gradPhi2Dphi2Dphi2D;
        hatG{2}(ind) = W * (gradPhi2Dphi2Dphi2D .* Q1);
        hatG{3}(ind) = W * (gradPhi2Dphi2Dphi2D .* Q2);
      end  % for i
    end  % for l
  end  % for j
end  % for m
end  % function