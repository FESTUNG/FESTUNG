nction [ hatM ] = computeHatM1D( p , ord )

global gPsi1D

[~, W] = gaussQuadRule1D(ord);
hatM = zeros(p+1, p+1);

for i = 1 : p+1
    for j = 1 : i
        hatM(i,j) = W * ( gPsi1D(:,i) .* gPsi1D(:,j) );
        hatM(j,i) = hatM(i,j);
    end  % for j
end  % for i

end  % function