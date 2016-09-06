function [ hatSoffdiag ] = computeHatSoffdiag( p , ord )

global gPhi1D

numBases = (p+1)^2;
[~, W] = gaussQuadRule1D(ord);
hatSoffdiag = zeros(numBases, numBases, 4);

for i = 1 : numBases
    for j = 1 : numBases
        hatSoffdiag(i,j,1) = W * ( gPhi1D(:,i,1) .* gPhi1D(:,j,2) );
        hatSoffdiag(i,j,2) = W * ( gPhi1D(:,i,2) .* gPhi1D(:,j,1) );
        hatSoffdiag(i,j,3) = W * ( gPhi1D(:,i,3) .* gPhi1D(:,j,4) );
        hatSoffdiag(i,j,4) = W * ( gPhi1D(:,i,4) .* gPhi1D(:,j,3) );
    end  % for j
end  % for i

end  % function