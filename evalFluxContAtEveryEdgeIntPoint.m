function dataEvalTest = evalFluxContAtEveryEdgeIntPoint(g, problemData, fluxCont, cCont, Nlambda)
validateattributes(fluxCont, {'function_handle'}, {}, mfilename, 'fluxCont');

p = Nlambda-1;  qOrd = max(2*p, 1);  [Q,~] = quadRule1D(qOrd);
K = g.numT;
Kedge = g.numE;

dataEvalTest = zeros( K, 3, 2, size(Q,2) );
F1 = @(X1, X2) g.B(:,1,1)*X1 + g.B(:,1,2)*X2 + g.coordV0T(:,1,1)*ones(size(X1));
F2 = @(X1, X2) g.B(:,2,1)*X1 + g.B(:,2,2)*X2 + g.coordV0T(:,1,2)*ones(size(X1));
for n = 1:3
    [Q1, Q2] = gammaMap(n, Q);
    Kkn = ~g.markE0Tint(:, n) .* g.markE0TbdrD(:, n);
    dataEvalTest(:, n, 1, :) = Kkn .* problemData.u1Cont( 0, F1(Q1, Q2), F2(Q1, Q2) ) .* cCont(:, :, n)  ;
    dataEvalTest(:, n, 2, :) = Kkn .* problemData.u2Cont( 0, F1(Q1, Q2), F2(Q1, Q2) ) .* cCont(:, :, n)  ;
end
%
lambdaEval = problemData.cDiscLambda(:,:) * problemData.basesOnGamma.phi1D{qOrd}';
for n = 1:3
    [Q1, Q2] = gammaMap(n, Q);
    Kkn = ~g.markE0TbdrD(:, n);
    dataEvalTest(:, n, 1, :) = problemData.u1Cont( 0, F1(Q1, Q2), F2(Q1, Q2) ) ...
        .* (sparse( 1:K, g.E0T(:, 1), Kkn .* ones(K, 1) , K, Kedge ) * lambdaEval(:,:)) ;
    dataEvalTest(:, n, 2, :) = problemData.u2Cont( 0, F1(Q1, Q2), F2(Q1, Q2) ) ...
        .* (sparse( 1:K, g.E0T(:, 1), Kkn .* ones(K, 1) , K, Kedge ) * lambdaEval(:,:)) ;
end


end % function