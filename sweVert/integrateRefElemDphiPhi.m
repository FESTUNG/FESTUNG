function [ hatHc, hatHx , hatHy ] = computeHatH( p , ord )

global gPhi2D gGradPhi2D

numBases = (p+1)^2;
[Q1, Q2, W] = gaussQuadRule2D(ord);
hatHc = zeros(numBases, numBases, 2);
hatHx = zeros(numBases, numBases, 2);
hatHy = zeros(numBases, numBases, 2);

for i = 1 : numBases
    for j = 1 : numBases
        for m = 1 : 2
            product = gPhi2D(:,j) .* gGradPhi2D(:,i,m);
            hatHc(i,j,m) = W * product;
            hatHx(i,j,m) = W * ( product .* Q1 );
            hatHy(i,j,m) = W * ( product .* Q2 );
        end  % for m
    end  % for j
end  % for i

end  % function