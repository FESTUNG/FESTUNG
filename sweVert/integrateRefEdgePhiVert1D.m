function [ hatIntV ] = computeHatIntV( p , ord )

global gPhi1D

numBases = p+1;
[~, W] = gaussQuadRule1D(ord);
hatIntV = zeros(numBases, numBases^2, 2);

for n = 1 : 2
    for j = 1 : numBases^2
        hatIntV(:,j,n) = (W * gPhi1D(:,j,n+2)) * ones(numBases,1);
    end  % for j
end  % for n

end  % function