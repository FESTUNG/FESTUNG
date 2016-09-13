function [ globJD ] = assembleGlobJDsub( g , p , ord , markE0Tbdr , cDAlg )

global gPhi1D

K = g.numTsub;
N = (p+1)^2;
[Q, W] = gaussQuadRule1D(ord);
Q2X = @(X,Y) g.deltaX * ones(g.numTsub,1) * X + g.coordV0Tsub(:,1,1) * ones(size(X));
Q2Y = @(X,Y) g.BAsub * X + g.DAsub * Y + g.ACBDsub * (X .* Y) + g.coordV0Tsub(:,1,2) * ones(size(X));

globJD = cell(2, 1);
globJD{1} = zeros(K, N);
globJD{2} = zeros(K, N);

for n = 1 : 4
    switch n
        case 1, QX = Q;              QY = zeros(size(Q));
        case 2, QX = Q;              QY = ones(size(Q));
        case 3, QX = ones(size(Q));  QY = Q;
        case 4, QX = zeros(size(Q)); QY = Q;
    end  % siwtch
    cDn = cDAlg( Q2X(QX.', QY.') , Q2Y(QX.', QY.') );
    for i = 1 : N
        integral = cDn * ( W .* gPhi1D(:,i,n)' )';
        globJD{1}(:,i) = globJD{1}(:,i) - integral .* markE0Tbdr(:,n) .* g.NuLengthE0Tsub(:,n,1);
        globJD{2}(:,i) = globJD{2}(:,i) - integral .* markE0Tbdr(:,n) .* g.NuLengthE0Tsub(:,n,2);
    end  % for i
end  % for n

globJD{1} = reshape(globJD{1}', K*N, 1);
globJD{2} = reshape(globJD{2}', K*N, 1);

end  % function