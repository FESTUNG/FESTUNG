function [ hatMc , hatMx ] = computeHatM( p , ord )

global gPhi2D

numBases = (p+1)^2;
[Q1, ~, W] = gaussQuadRule2D(ord);
hatMc = zeros(numBases,numBases);
hatMx = zeros(numBases,numBases);

for i = 1 : numBases
    for j = 1 : i
        product = gPhi2D(:,i) .* gPhi2D(:,j);
        hatMc(i,j) = W * product;
        hatMc(j,i) = hatMc(i,j);
        hatMx(i,j) = W * ( product .* Q1 );
        hatMx(j,i) = hatMx(i,j);
    end  % for j
end  % for i

end  % function