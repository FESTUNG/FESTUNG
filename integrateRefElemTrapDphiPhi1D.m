function ret = integrateRefElemTrapDphiPhi1D(N, qOrd, basesOnQuad2D, basesOnQuad1D)
validateattributes(N, {'numeric'}, {'numel', 2}, mfilename, 'N')
validateattributes(basesOnQuad2D, {'struct'}, {}, mfilename, 'basesOnQuad2D')
validateattributes(basesOnQuad1D, {'struct'}, {}, mfilename, 'basesOnQuad1D')
[Q, W] = quadRule1D(qOrd); [Q1, Q2] = meshgrid(Q); W = W' * W;
Q1 = Q1(:); Q2 = Q2(:); W = W(:)';
ret = { zeros(N(1), N(2), 2), zeros(N(1), N(2), 2), zeros(N(1), N(2), 2) };
for m = 1 : 2
  for i = 1 : N(1)
    for j = 1 : N(2)
      gradPhi2Dphi1D = basesOnQuad2D.gradPhi2D(:,i,m) .* kron(basesOnQuad1D.phi1D(:,j), ones(size(Q)).');
      ret{1}(i,j,m) = W * gradPhi2Dphi1D;
      ret{2}(i,j,m) = W * ( gradPhi2Dphi1D .* Q1 );
      ret{3}(i,j,m) = W * ( gradPhi2Dphi1D .* Q2 );
    end % for j
  end % for i
end % for m
end % function
