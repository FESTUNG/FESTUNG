%
function ret = assembleVecEdgeMuFuncContVal(g, markE0Tbdr, funcCont, Nmu, basesOnGamma )

KEdge = g.numE;

% Determine quadrature rule
p = Nmu - 1;
qOrd = 2*p+1;
[Q, W] = quadRule1D(qOrd);
R = size( Q, 2 );
ret = zeros(KEdge, Nmu);
Q2X1 = @(X1,X2) g.B(:,1,1)*X1 + g.B(:,1,2)*X2 + g.coordV0T(:,1,1)*ones(size(X1));
Q2X2 = @(X1,X2) g.B(:,2,1)*X1 + g.B(:,2,2)*X2 + g.coordV0T(:,1,2)*ones(size(X1));
for n = 1 : 3
    [Q1, Q2] = gammaMap(n, Q);
    funcOnQuad = funcCont(Q2X1(Q1, Q2), Q2X2(Q1, Q2));
    Kkn = markE0Tbdr(:, n) .* g.areaE0T(:,n);
    %Matlab 2016b
%     ret(g.E0T(:, n),:) = ret(g.E0T(:, n),:) + Kkn .* funcOnQuad * ( W' .* basesOnGamma.phi1D{qOrd}(:,:) );
    %Matlab 2016a
    ret(g.E0T(:, n),:) = ret(g.E0T(:, n),:) + repmat(Kkn, 1, R) .* funcOnQuad * ( repmat(W', 1, Nmu) .* basesOnGamma.phi1D{qOrd}(:,:) );
end
ret = reshape(ret',KEdge*Nmu,1);
end

