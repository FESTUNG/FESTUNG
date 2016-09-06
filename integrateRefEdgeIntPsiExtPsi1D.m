function [ hatQoffdiag1D ] = computeHatQoffdiag1D( p )

global gPsi1Dbnd

numBases = p+1;
hatQoffdiag1D = zeros(numBases, numBases, 2);

for i = 1 : numBases
    for j = 1 : numBases
        hatQoffdiag1D(i,j,1) = gPsi1Dbnd(1,i) * gPsi1Dbnd(2,j);
        hatQoffdiag1D(i,j,2) = gPsi1Dbnd(2,i) * gPsi1Dbnd(1,j);
    end  % for j
end  % for i

end  % function