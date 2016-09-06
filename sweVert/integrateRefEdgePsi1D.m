function [ hatR1D ] = computeHatR1D( p )

global gPsi1Dbnd

numBases = p+1;
hatR1D = zeros(numBases, numBases, 2);

for j = 1 : numBases
    for n = 1 : 2
        hatR1D(:,j,n) = gPsi1Dbnd(n,j) * ones(numBases,1);
    end  % for n
end  % for j

end  % function