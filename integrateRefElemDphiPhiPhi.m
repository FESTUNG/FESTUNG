function [ hatGc , hatGx , hatGy ] = computeHatG( p , ord )

global gPhi2D gGradPhi2D

numBases = (p+1)^2;
[Q1, Q2, W] = gaussQuadRule2D(ord);
hatGc = zeros(numBases, numBases, numBases, 2);
hatGx = zeros(numBases, numBases, numBases, 2);
hatGy = zeros(numBases, numBases, numBases, 2);

for i = 1 : numBases
    for j = 1 : numBases
        for k = 1 : j
            for m = 1 : 2
                product = gGradPhi2D(:,i,m) .* gPhi2D(:,j) .* gPhi2D(:,k);
                hatGc(i,j,k,m) = W * product;
                hatGc(i,k,j,m) = hatGc(i,j,k,m);
                hatGx(i,j,k,m) = W * (product .* Q1);
                hatGx(i,k,j,m) = hatGx(i,j,k,m);
                hatGy(i,j,k,m) = W * (product .* Q2);
                hatGy(i,k,j,m) = hatGy(i,j,k,m);
            end  % for m
        end  % for k
    end  % for j
end  % for i

end  % function