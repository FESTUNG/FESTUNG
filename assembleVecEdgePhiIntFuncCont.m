function [ globKD ] = assembleGlobKDsub( g , p , ord , markE0Tbdr , cDAlg , eta )

global gPhi1D

K = g.numTsub;
N = (p+1)^2;
[Q, W] = gaussQuadRule1D(ord);
Q2X = @(X,Y) g.deltaX * ones(g.numTsub,1) * X + g.coordV0Tsub(:,1,1) * ones(size(X));
Q2Y = @(X,Y) g.BAsub * X + g.DAsub * Y + g.ACBDsub * (X .* Y) + g.coordV0Tsub(:,1,2) * ones(size(X));

globKD = zeros(K, N);

for n = 1 : 4
    switch n
        case 1, QX = Q;              QY = zeros(size(Q));
        case 2, QX = Q;              QY = ones(size(Q));
        case 3, QX = ones(size(Q));  QY = Q;
        case 4, QX = zeros(size(Q)); QY = Q;
    end  % siwtch
    cDn = cDAlg( Q2X(QX.', QY.') , Q2Y(QX.', QY.') );
    for i = 1 : N
        globKD(:,i) = globKD(:,i) + eta *  markE0Tbdr(:,n) .* ( cDn * ( W' .* gPhi1D(:,i,n) ) );
    end  % for i
end  % for n

globKD = reshape(globKD', K*N, 1);

end  % function