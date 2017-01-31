function basesOnQuad = computeBasesOnQuad1D(p, qOrd, basesOnQuad)
[Q, ~] = quadRule1D(qOrd);
R = length(Q); N = (p+1);
basesOnQuad.phi0D = zeros(N, 2);
basesOnQuad.phi1D = zeros(R, N);
basesOnQuad.gradPhi1D = zeros(R, N);
for i = 1 : N
  basesOnQuad.phi1D(:, i) = phi1D(i, Q);
  basesOnQuad.gradPhi1D(:, i) = gradPhi1D(i, Q);
  for n = 1 : 2
    basesOnQuad.phi0D(i, n) = phi1D(i, n-1);
  end % for n
end % for i
end  % function