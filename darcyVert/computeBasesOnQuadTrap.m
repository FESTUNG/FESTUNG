function basesOnQuad = computeBasesOnQuadTrap(p, qOrd, basesOnQuad)
[Q, ~] = quadRule1D(qOrd);
[Q1, Q2] = meshgrid(Q); Q1 = Q1(:)'; Q2 = Q2(:)';
R1D = length(Q); R2D = length(Q1); N = (p+1) * (p+1);

basesOnQuad.phi1D = zeros(R1D, N, 4);
basesOnQuad.phi2D = zeros(R2D, N);
basesOnQuad.gradPhi2D = zeros(R2D, N, 2);

for i = 1 : N
  basesOnQuad.phi2D(:,i) = phiTrap(i, Q1, Q2);
  for m = 1 : 2
    basesOnQuad.gradPhi2D(:,i,m) = gradPhiTrap(i, m, Q1, Q2);
  end % for m
  for n = 1 : 2
    gamma_n = n - 1;
    basesOnQuad.phi1D(:,i,n) = phiTrap(i, Q, gamma_n);  % edges 1 (bottom), 2 (top)
    basesOnQuad.phi1D(:,i,5-n) = phiTrap(i, gamma_n, Q); % edges 3 (right), 4 (left)
  end % for n
end  % for i

% basesOnQuad.psi1D = zeros(length(Q), p);
% basesOnQuad.gDerivPsi1D = zeros(length(Q), p);
% basesOnQuad.psi1Dbnd = zeros(2, p);
% 
% for i = 1 : p+1
%     basesOnQuad.psi1D(:,i) = psi1D(i, Q);
%     basesOnQuad.gDerivPsi1D(:,i) = derivPsi(i, Q);
%     basesOnQuad.psi1Dbnd(1,i) = psi1D(i,1);
%     basesOnQuad.psi1Dbnd(2,i) = psi1D(i,0);
% end  % for

end  % function