function [ PhiPhioffdiag ] = computePhiPhioffdiag( p , ord )

global gPhi1D

numBases = (p+1)^2;
[~, W] = gaussQuadRule1D(ord); 
PhiPhioffdiag = zeros(numBases, numBases, length(W), 2);


for i = 1 : numBases
    for j = 1 : numBases
        PhiPhioffdiag(i, j, :, 1) = gPhi1D(:,i,3) .* gPhi1D(:,j,4) .* W';
        PhiPhioffdiag(i, j, :, 2) = gPhi1D(:,i,4) .* gPhi1D(:,j,3) .* W';
    end  % for i
end  % for j

end  % function