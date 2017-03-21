function ret = integrateRefElemTrapDphiPhi(N, qOrd, basesOnQuad)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
[Q, W] = quadRule1D(qOrd); [Q1, Q2] = meshgrid(Q); W = W' * W;
Q1 = Q1(:); Q2 = Q2(:); W = W(:)';
ret = { zeros(N, N, 2), zeros(N, N, 2), zeros(N, N, 2) };
for m = 1 : 2
  for j = 1 : N
    for i = 1 : N
      phi2DgradPhi2D = basesOnQuad.phi2D{qOrd}(:,j) .* basesOnQuad.gradPhi2D{qOrd}(:,i,m);
      ret{1}(i,j,m) = W * phi2DgradPhi2D;
      ret{2}(i,j,m) = W * ( phi2DgradPhi2D .* Q1 );
      ret{3}(i,j,m) = W * ( phi2DgradPhi2D .* Q2 );
    end  % for i
  end  % for j
end  % for m
end  % function