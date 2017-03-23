function dataEval = evalUContAtEveryEdgeIntPoint(g, fCont, Nmu)
validateattributes(fCont, {'function_handle'}, {}, mfilename, 'fCont');

p = Nmu-1;  qOrd = max(2*p+1, 1);  [Q1,~] = quadRule1D(qOrd);
K = g.numT;

dataEval = cell(2,1);
dataEval{1} = zeros( size(Q1,2), K, 3 );
dataEval{2} = zeros( size(Q1,2), K, 3 );

F1 = @(X1, X2) g.B(:,1,1)*X1 + g.B(:,1,2)*X2 + g.coordV0T(:,1,1)*ones(size(X1));
F2 = @(X1, X2) g.B(:,2,1)*X1 + g.B(:,2,2)*X2 + g.coordV0T(:,1,2)*ones(size(X1));
for n = 1:3
    [x1, x2] = gammaMap( n, Q1 );
    [dataEval{1}(:, :, n), dataEval{2}(:, :, n)]= fCont( F1(x1, x2)', F2(x1, x2)' );
end

end % function