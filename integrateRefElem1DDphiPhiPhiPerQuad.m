function ret = integrateRefElem1DDphiPhiPhi(N, qOrd, basesOnQuad)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
[Q, W] = quadRule1D(qOrd);
ret = { zeros(N), zeros(N) }; % [N x N]
for i = 1 : N
  for j = 1 : N
    for l = 1 : N
      ret{1}(i, j, l) = W * ( basesOnQuad.gradPhi1D(:, i) .* basesOnQuad.phi1D(:, j) .* basesOnQuad.phi1D(:, l) );
      ret{2}(i, j, l) = (W .* Q) * ( basesOnQuad.gradPhi1D(:, i) .* basesOnQuad.phi1D(:, j) .* basesOnQuad.phi1D(:, l) );
    end % for l
  end % for j
end % for i
end % function