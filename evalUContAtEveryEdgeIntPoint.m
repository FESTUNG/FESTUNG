function dataEval = evalUContAtEveryEdgeIntPoint(g, fluxCont, Nlambda)
validateattributes(fluxCont, {'function_handle'}, {}, mfilename, 'fluxCont');

p = Nlambda-1;  qOrd = max(2*p, 1);  [Q1,~] = quadRule1D(qOrd);
Kedge = g.numE;

dataEval = zeros( 2, size(Q1,2), Kedge );

for iE = 1:Kedge
    localIdx = g.E0E(iE,1);
    adjTri = g.T0E(iE, 1);
    
    [x1, x2] = gammaMap( localIdx, Q1 );
    
    F1 = @(X1, X2) g.B(adjTri,1,1)*X1 + g.B(adjTri,1,2)*X2 + g.coordV0T(adjTri,1,1)*ones(size(X1));
    F2 = @(X1, X2) g.B(adjTri,2,1)*X1 + g.B(adjTri,2,2)*X2 + g.coordV0T(adjTri,1,2)*ones(size(X1));

    dataEval(:, :, iE)  = fluxCont( F1(x1, x2), F2(x1, x2) );
end
end % function