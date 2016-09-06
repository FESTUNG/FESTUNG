function [ globAlpha ] = assembleGlobAlpha( g , p , ord , uDG , hDG , gNum , ...
    PhiPhiDiag , PhiPhiOffdiag , markE0Tint)

global gPsi1Dbnd gPhi1D

[K, ~] = size(uDG);

[~, W] = gaussQuadRule1D(ord);
R = length(W);
markE0Tint3 = markE0Tint(1:end-1,3);

globAlpha = sparse( g.numTsupra * (p+1)^2 , g.numTsupra * (p+1)^2 );

uLeft   = zeros(g.numTsupra * R , 1);
uRight  = zeros(g.numTsupra * R , 1);
hLeft   = zeros(g.NX , 1);
hRight  = zeros(g.NX , 1);

for i = 1 : p+1
    hLeft   = hLeft     + kron(hDG(:,i), gPsi1Dbnd(1,i));
    hRight  = hRight    + kron(hDG(:,i), gPsi1Dbnd(2,i));
end  % for

for i = 1 : (p+1)^2
    uLeft   = uLeft     + kron(uDG(:,i), gPhi1D(:,i,3));
    uRight  = uRight    + kron(uDG(:,i), gPhi1D(:,i,4));
end  % for

uLeft = uLeft(1:(g.numTsupra-1)*R);
uRight = uRight(R+1:(g.numTsupra)*R);

hatU = 0.5 * (uLeft + uRight);

hLeft = hLeft(1:g.NX-1);
hRight = hRight(2:g.NX);

hatH = kron( 0.5 * (hLeft + hRight) , ones(R,1) );
hatH = [hatH; zeros(R,1)];
hatH = repmat(hatH, g.NZsupra, 1);
hatH = hatH(1:length(hatU));

% Check: hatH should not contain negative values
hatH(hatH < 10^(-6)) = 10^(-6);

lambdaMax = 0.75 * abs(hatU) + 0.25 * sqrt( hatU.*hatU + 4 * gNum * hatH );
lambdaMax = reshape(lambdaMax, R, g.numTsupra-1).';
lambdaMax = bsxfun(@times, lambdaMax, markE0Tint3);

for r = 1 : R
    length3 = g.lengthE0Tsupra(:,3) .* [lambdaMax(:,r); 0];
    length4 = g.lengthE0Tsupra(:,4) .* [0; lambdaMax(:,r)];
    
    assert(isreal(length3), 'Complex values occur in length3.');
    assert(isreal(length4), 'Complex values occur in length4.');

    globAlpha = globAlpha + kron( spdiags(length3, 0, K, K), PhiPhiDiag(:,:,r,1) ) ...
        + kron( spdiags(length4, 0, K, K), PhiPhiDiag(:,:,r,2) ) ...
        - kron( bsxfun(@times, g.markE0TE0Tsupra{3}, length3), PhiPhiOffdiag(:,:,r,1) ) ...
        - kron( bsxfun(@times, g.markE0TE0Tsupra{4}, length4), PhiPhiOffdiag(:,:,r,2) );
    
    assert(isreal(globAlpha), 'Complex values occur in globAlpha.');
    
end  % for

end  % function