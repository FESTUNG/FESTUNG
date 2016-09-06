function [ hatSoffdiag ] = integrateRefEdgePhiIntPhiExt( p , ord , basesOnQuad )

numBases = (p+1)^2;
[~, W] = gaussQuadRule1D(ord);
hatSoffdiag = zeros(numBases, numBases, 4);

for i = 1 : numBases
    for j = 1 : numBases
        hatSoffdiag(i,j,1) = W * ( basesOnQuad.gPhi1D(:,i,1) .* basesOnQuad.gPhi1D(:,j,2) );
        hatSoffdiag(i,j,2) = W * ( basesOnQuad.gPhi1D(:,i,2) .* basesOnQuad.gPhi1D(:,j,1) );
        hatSoffdiag(i,j,3) = W * ( basesOnQuad.gPhi1D(:,i,3) .* basesOnQuad.gPhi1D(:,j,4) );
        hatSoffdiag(i,j,4) = W * ( basesOnQuad.gPhi1D(:,i,4) .* basesOnQuad.gPhi1D(:,j,3) );
    end  % for j
end  % for i

end  % function