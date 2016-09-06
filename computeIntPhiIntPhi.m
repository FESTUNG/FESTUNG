function [ PhiPhidiag ] = computePhiPhidiag( p , ord )

global gPhi1D

numBases = (p+1)^2;
[~, W] = gaussQuadRule1D(ord);
PhiPhidiag = zeros(numBases, numBases, length(W), 2);

for n = 1 : 2
  for i = 1 : numBases
    for j = 1 : i
      PhiPhidiag(i, j, :, n) = gPhi1D(:, i, n+2) .* gPhi1D(:, j, n+2) .* W';
      PhiPhidiag(j, i, :, n) = PhiPhidiag(i,j,:,n);
    end  % for
  end  % for
end  % for

end  % function