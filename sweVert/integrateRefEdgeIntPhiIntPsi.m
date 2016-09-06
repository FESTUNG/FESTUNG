function [ hatCdiag ] = computeHatCdiag( p , ord )

global gPsi1Dbnd gPsi1D gPhi1D

numBases = p+1;
[~, W] = gaussQuadRule1D(ord);
hatCdiag = zeros(numBases^2, numBases, 4);

for n = 1 : 2
    for i = 1 : numBases^2
        for j = 1 : numBases
            hatCdiag(i,j,n) = W * ( gPhi1D(:,i,n) .* gPsi1D(:,j) );
        end  % for j
    end  % for i
end  % for n

for n = 3 : 4
    for i = 1 : numBases^2
        for j = 1 : numBases
            hatCdiag(i,j,n) = (W * gPhi1D(:,i,n)) * gPsi1Dbnd(n-2,j);
        end  % for j
    end  % for i
end  % for n

end  % function