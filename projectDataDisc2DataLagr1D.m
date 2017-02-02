function dataLagr = projectDataDisc2DataLagr1D(dataDisc, p)
[K, N] = size(dataDisc);
if nargin == 1
  p = (sqrt(8*N+1)-3)/2;
end % if
switch p
  case 0,    L1 = 1/2;       % locally constant
  case 1,    L1 = [0, 1];           % locally linear
  otherwise, L1 = [0, 1/2, 1]; % locally quadratic
end % switch
dataLagr = zeros(K, length(L1));
for i = 1 : N
  dataLagr = dataLagr + dataDisc(:, i) * phi1D(i, L1);
end % for
end % function