function [ hatRdiag ] = computeHatRdiag( p , ord )

global gPhi1D

numBases = (p+1)^2;
[~, W] = gaussQuadRule1D(ord);
hatRdiag = zeros(numBases, numBases, numBases, 4);

for n = 1 : 4  % number of edge
    for s = 1 : numBases  % basis-function of included function D(t)
        for i = 1 : s
            for j = 1 : i
                hatRdiag(i,j,s,n) = W * ( gPhi1D(:,i,n) .* gPhi1D(:,s,n) .* gPhi1D(:,j,n) );
                hatRdiag(j,i,s,n) = hatRdiag(i,j,s,n);
                hatRdiag(i,s,j,n) = hatRdiag(i,j,s,n);
                hatRdiag(j,s,i,n) = hatRdiag(i,j,s,n);
                hatRdiag(s,i,j,n) = hatRdiag(i,j,s,n);
                hatRdiag(s,j,i,n) = hatRdiag(i,j,s,n);
            end  % for j
        end  % for i
    end  % for s
end  % for n

end  % function