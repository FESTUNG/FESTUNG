function [ hatQdiag1D ] = computeHatQdiag1D( p )

global gPsi1Dbnd

numBases = p+1;
hatQdiag1D = zeros(numBases, numBases, 2);

for n = 1 : 2
    for i = 1 : numBases
        for j = 1 : i
            hatQdiag1D(i,j,n) = gPsi1Dbnd(n,i) * gPsi1Dbnd(n,j);
            hatQdiag1D(j,i,n) = hatQdiag1D(i,j,n);
        end  % for j
    end  % for i
end  % for n

end  % function