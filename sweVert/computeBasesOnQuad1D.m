function basesOnQuad = computeBasesOnQuad1D(p, basesOnQuad, requiredOrders)
if nargin < 3
  if p > 0
    requiredOrders = [2*p, 2*p+1]; 
  else
    requiredOrders = 1; 
  end % if
end % if
%
N = p+1;

basesOnQuad.phi0D = cell(max(requiredOrders), 1);
basesOnQuad.phi1D = cell(max(requiredOrders), 1);
basesOnQuad.gradPhi1D = cell(max(requiredOrders), 1);

for qOrd = requiredOrders
  [Q, ~] = quadRule1D(qOrd);  R = length(Q); 
  basesOnQuad.phi0D{qOrd} = zeros(N, 2);
  basesOnQuad.phi1D{qOrd} = zeros(R, N);
  basesOnQuad.gradPhi1D{qOrd} = zeros(R, N);
  for i = 1 : N
    basesOnQuad.phi1D{qOrd}(:, i) = phi1D(i, Q);
    basesOnQuad.gradPhi1D{qOrd}(:, i) = gradPhi1D(i, Q);
    for n = 1 : 2
      basesOnQuad.phi0D{qOrd}(i, n) = phi1D(i, n-1);
    end % for n
  end % for i
end % for qOrd
end  % function