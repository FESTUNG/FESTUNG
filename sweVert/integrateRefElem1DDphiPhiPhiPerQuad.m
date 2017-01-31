function ret = integrateRefElem1DDphiPhiPhiPerQuad(N, qOrd, basesOnQuad)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
[Q, W] = quadRule1D(qOrd); R = length(W);
ret = { zeros(N, N, N, R), zeros(N, N, N, R), zeros(N, N, N, R) }; % [N x N x N x R]
for i = 1 : N
  for j = 1 : N
    for l = 1 : N
      ret{1}(i, j, l, :) = W.' .* ( basesOnQuad.gradPhi1D(:, i) .* basesOnQuad.phi1D(:, j) .* basesOnQuad.phi1D(:, l) );
      ret{2}(i, j, l, :) = (W .* Q).' .* ( basesOnQuad.gradPhi1D(:, i) .* basesOnQuad.phi1D(:, j) .* basesOnQuad.phi1D(:, l) );
    end % for l
  end % for j
end % for i
end % function