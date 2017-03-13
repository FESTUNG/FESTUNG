function dataEval = evalSourceContAtEveryIntPoint(g, sourceCont, N)
validateattributes(sourceCont, {'function_handle'}, {}, mfilename, 'sourceCont');
K = g.numT;
p = (sqrt(8*N+1)-3)/2;  qOrd = max(2*p,1);  [Q1, Q2, ~] = quadRule2D(qOrd);

F1 = @(X1, X2) g.B(:,1,1)*X1 + g.B(:,1,2)*X2 + g.coordV0T(:,1,1)*ones(size(X1));
F2 = @(X1, X2) g.B(:,2,1)*X1 + g.B(:,2,2)*X2 + g.coordV0T(:,1,2)*ones(size(X1));

% dataEval = zeros( K, size(Q1, 2) );
dataEval = sourceCont(F1(Q1, Q2), F2(Q1, Q2));
end % function