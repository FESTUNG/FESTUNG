function [ gradPsiPhiC , gradPsiPhiX , gradPsiPhiY ] = computeGradPsiPhi( p , ord )

global gPhi2D gDerivPsi1D

numBases = p+1;
[Q1, Q2, W] = gaussQuadRule2D(ord);
helper = length(gPhi2D(:,1)) / length(gDerivPsi1D(:,1));

gradPsiPhiC = zeros(numBases, numBases^2, length(W));
gradPsiPhiX = zeros(numBases, numBases^2, length(W));
gradPsiPhiY = zeros(numBases, numBases^2, length(W));

for i = 1 : numBases
    for j = 1 : numBases^2
        product = kron(gDerivPsi1D(:,i), ones(helper,1)) .* gPhi2D(:,j) .* W';
        gradPsiPhiC(i,j,:) = product;
        gradPsiPhiX(i,j,:) = product .* Q1;
        gradPsiPhiY(i,j,:) = product .* Q2;
    end  % for j
end  % for i

end  % function