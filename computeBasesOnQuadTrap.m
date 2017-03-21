function basesOnQuad = computeBasesOnQuadTrap(p, basesOnQuad, requiredOrders)
if nargin < 3
  if p > 0
    requiredOrders = [2*p, 2*p+1]; 
  else
    requiredOrders = 1; 
  end % if
end % if
%
N = (p+1) * (p+1);

basesOnQuad.phi1D = cell(max(requiredOrders), 1);
basesOnQuad.phi2D = cell(max(requiredOrders), 1);
basesOnQuad.gradPhi2D = cell(max(requiredOrders), 1);

for qOrd = requiredOrders
  [Q, ~] = quadRule1D(qOrd);
  [Q1, Q2] = meshgrid(Q); Q1 = Q1(:)'; Q2 = Q2(:)';
  R1D = length(Q); R2D = length(Q1);

  basesOnQuad.phi1D{qOrd} = zeros(R1D, N, 4);
  basesOnQuad.phi2D{qOrd} = zeros(R2D, N);
  basesOnQuad.gradPhi2D{qOrd} = zeros(R2D, N, 2);
  for i = 1 : N
    basesOnQuad.phi2D{qOrd}(:,i) = phiTensorProduct(i, Q1, Q2, @phi1D, @phi1D);
    for m = 1 : 2
      basesOnQuad.gradPhi2D{qOrd}(:,i,m) = gradPhiTensorProduct(i, m, Q1, Q2, @phi1D, @phi1D, @gradPhi1D, @gradPhi1D);
    end % for m
    for n = 1 : 2
      gamma_n = n - 1;
      basesOnQuad.phi1D{qOrd}(:,i,n) = phiTensorProduct(i, Q, gamma_n, @phi1D, @phi1D);  % edges 1 (bottom), 2 (top)
      basesOnQuad.phi1D{qOrd}(:,i,5-n) = phiTensorProduct(i, gamma_n, Q, @phi1D, @phi1D); % edges 3 (right), 4 (left)
    end % for n
  end  % for i
end % for qOrd
end  % function