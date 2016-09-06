function [ hatIntW ] = computeHatIntW( p , ord )

global gPhi1D

numBases = p+1;
[~, W] = gaussQuadRule1D(ord);
hatIntW = zeros(numBases, numBases^2, numBases^2, 2);

for n = 1 : 2
    for j = 1 : numBases^2
        for s = 1 : j
            hatIntW(:,j,s,n) = ( W * (gPhi1D(:,j,n+2) .* gPhi1D(:,s,n+2)) ) * ones(numBases,1);
            hatIntW(:,s,j,n) = hatIntW(:,j,s,n);
        end  % for s
    end  % for j
end  % for n

end  % function