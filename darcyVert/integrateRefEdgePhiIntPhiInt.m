function [ hatSdiag ] = integrateRefEdgePhiIntPhiInt( p , ord , basesOnQuad )

numBases = (p+1)^2;
[~, W] = gaussQuadRule1D(ord);
hatSdiag = zeros(numBases, numBases, 4);

for n = 1 : 4  % number edge
    for i = 1 : numBases
        for j = 1 : i
            hatSdiag(i,j,n) = W * ( basesOnQuad.gPhi1D(:,i,n) .* basesOnQuad.gPhi1D(:,j,n) );
            hatSdiag(j,i,n) = hatSdiag(i,j,n);
        end  % for j
    end  % for i
end  % for n

end  % function