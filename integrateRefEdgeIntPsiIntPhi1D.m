function [ hatIntQdiag ] = computeHatIntQdiag( p , ord )

global gPhi1D gPsi1Dbnd

numBases = p+1;
[~, W] = gaussQuadRule1D(ord);
hatIntQdiag = zeros(numBases, numBases^2, 2);


for n = 1 : 2
    for i = 1 : numBases
        for j = 1 : numBases^2
            hatIntQdiag(i,j,n) = W * ( gPsi1Dbnd(n,i) .* gPhi1D(:,j,n+2) );
        end  % for j
    end  % for i
end  % for n

end  % function