function [ hatRoffdiag ] = computeHatRoffdiag( p , ord )

global gPhi1D

numBases = (p+1)^2;
[~, W] = gaussQuadRule1D(ord);
hatRoffdiag = zeros(numBases, numBases, numBases, 4);

for s = 1 : numBases  % basis-function of included function D(t)
    for i = 1 : numBases
        for j = 1 : s
            hatRoffdiag(i,j,s,1) = W * ( gPhi1D(:,i,1) .* gPhi1D(:,s,2) .* gPhi1D(:,j,2) );
            hatRoffdiag(i,s,j,1) = hatRoffdiag(i,j,s,1);
            hatRoffdiag(i,j,s,2) = W * ( gPhi1D(:,i,2) .* gPhi1D(:,s,1) .* gPhi1D(:,j,1) );
            hatRoffdiag(i,s,j,2) = hatRoffdiag(i,j,s,2);
            hatRoffdiag(i,j,s,3) = W * ( gPhi1D(:,i,3) .* gPhi1D(:,s,4) .* gPhi1D(:,j,4) );
            hatRoffdiag(i,s,j,3) = hatRoffdiag(i,j,s,3);
            hatRoffdiag(i,j,s,4) = W * ( gPhi1D(:,i,4) .* gPhi1D(:,s,3) .* gPhi1D(:,j,3) );
            hatRoffdiag(i,s,j,4) = hatRoffdiag(i,j,s,4);
        end  % for j
    end  % for i
end  % for s

end  % function