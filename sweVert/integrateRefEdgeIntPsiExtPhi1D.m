function [ hatIntQoffdiag ] = computeHatIntQoffdiag( p , ord )

global gPhi1D gPsi1Dbnd

numBases = p+1;
[~, W] = gaussQuadRule1D(ord);
hatIntQoffdiag = zeros(numBases, numBases^2, 2);

for i = 1 : numBases
    for j = 1 : numBases^2
        hatIntQoffdiag(i,j,1) = W * ( gPsi1Dbnd(1,i) .* gPhi1D(:,j,4) );
        hatIntQoffdiag(i,j,2) = W * ( gPsi1Dbnd(2,i) .* gPhi1D(:,j,3) );
    end  % for j
end  % for i

end  % function