function [ hatEc , hatEx , hatEy ] = computeHatE( p , ord )

global gGradPhi2D gPsi1D

numBases = p+1;
[Q1, Q2, W] = gaussQuadRule2D(ord);
hatEc = zeros(numBases^2, numBases, 2);
hatEx = zeros(numBases^2, numBases, 2);
hatEy = zeros(numBases^2, numBases, 2);

helper = length(gGradPhi2D(:,1,1)) / length(gPsi1D(:,1));

for m = 1 : 2
    for i = 1 : numBases^2
        for j = 1 : numBases
            product = kron(gPsi1D(:,j), ones(helper,1)) .* gGradPhi2D(:,i,m);
            hatEc(i,j,m) = W * product;
            hatEx(i,j,m) = W * (product .* Q1);
            hatEy(i,j,m) = W * (product .* Q2);
        end  % for j
    end  % for i
end  % for m

end  % function