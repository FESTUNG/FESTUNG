function dataEval = evalFuncContAtEveryEdgeIntPoint(g, funcCont, Nmu)
validateattributes(funcCont, {'function_handle'}, {}, mfilename, 'funcCont');
p = Nmu-1;  qOrd = max(2*p, 1);  [Q1,~] = quadRule1D(qOrd);
K = g.numT;
dataEval = zeros( K, size(Q1,2), 3 );

F1 = @(X1, X2) g.B(:,1,1)*X1 + g.B(:,1,2)*X2 + g.coordV0T(:,1,1)*ones(size(X1));
F2 = @(X1, X2) g.B(:,2,1)*X1 + g.B(:,2,2)*X2 + g.coordV0T(:,1,2)*ones(size(X1));
for n = 1:3
    [x1, x2] = gammaMap( n, Q1 );
    dataEval(:, :, n)  = funcCont( F1(x1, x2), F2(x1, x2) ) ;
end
end % function