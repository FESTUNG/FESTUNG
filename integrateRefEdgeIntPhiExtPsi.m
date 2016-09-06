function [ hatCoffdiag ] = computeHatCoffdiag( p , ord )

global gPsi1D gPsi1Dbnd gPhi1D

numBases = p+1;
[~, W] = gaussQuadRule1D(ord);
hatCoffdiag = zeros(numBases^2, numBases, 4);

for i = 1 : numBases^2
    for j = 1 : numBases
        hatCoffdiag(i,j,1) = W * ( gPhi1D(:,i,1) .* gPsi1D(:,j) );
        hatCoffdiag(i,j,2) = W * ( gPhi1D(:,i,2) .* gPsi1D(:,j) );
        hatCoffdiag(i,j,3) = (W * gPhi1D(:,i,3)) * gPsi1Dbnd(2,j);
        hatCoffdiag(i,j,4) = (W * gPhi1D(:,i,4)) * gPsi1Dbnd(1,j);
    end  % for j
end  % for i

end  % function