function [ hatRoffdiag ] = integrateRefEdgePhiIntPhiExtPhiExt( p , ord , basesOnQuad )

numBases = (p+1)^2;
[~, W] = gaussQuadRule1D(ord);
hatRoffdiag = zeros(numBases, numBases, numBases, 4);

for s = 1 : numBases  % basis-function of included function D(t)
    for i = 1 : numBases
        for j = 1 : s
            hatRoffdiag(i,j,s,1) = W * ( basesOnQuad.gPhi1D(:,i,1) .* basesOnQuad.gPhi1D(:,s,2) .* basesOnQuad.gPhi1D(:,j,2) );
            hatRoffdiag(i,s,j,1) = hatRoffdiag(i,j,s,1);
            hatRoffdiag(i,j,s,2) = W * ( basesOnQuad.gPhi1D(:,i,2) .* basesOnQuad.gPhi1D(:,s,1) .* basesOnQuad.gPhi1D(:,j,1) );
            hatRoffdiag(i,s,j,2) = hatRoffdiag(i,j,s,2);
            hatRoffdiag(i,j,s,3) = W * ( basesOnQuad.gPhi1D(:,i,3) .* basesOnQuad.gPhi1D(:,s,4) .* basesOnQuad.gPhi1D(:,j,4) );
            hatRoffdiag(i,s,j,3) = hatRoffdiag(i,j,s,3);
            hatRoffdiag(i,j,s,4) = W * ( basesOnQuad.gPhi1D(:,i,4) .* basesOnQuad.gPhi1D(:,s,3) .* basesOnQuad.gPhi1D(:,j,3) );
            hatRoffdiag(i,s,j,4) = hatRoffdiag(i,j,s,4);
        end  % for j
    end  % for i
end  % for s

end  % function