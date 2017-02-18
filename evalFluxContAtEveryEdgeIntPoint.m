function ret = evalFluxContAtEveryEdgeIntPoint(g, markE0T, fluxCont, cCont, Nlambda)
validateattributes(fluxCont, {'function_handle'}, {}, mfilename, 'fluxCont');

p = Nlambda-1;  qOrd = max(2*p, 1);  [Q,~] = quadRule1D(qOrd);
K = g.numT;

ret = zeros( K, 3, 2, size(Q,2) );
F1 = @(X1, X2) g.B(:,1,1)*X1 + g.B(:,1,2)*X2 + g.coordV0T(:,1,1)*ones(size(X1));
F2 = @(X1, X2) g.B(:,2,1)*X1 + g.B(:,2,2)*X2 + g.coordV0T(:,1,2)*ones(size(X1));
for n = 1:3
    [Q1, Q2] = gammaMap(n, Q);
    Kkn = markE0T(:, n);
 
% Das habe ich zuletzt verwendet
    [ret(:, n, 1, :), ret(:, n, 2, :)] = fluxCont( F1(Q1, Q2), F2(Q1, Q2), cCont(:, :, n) );
    ret(:, n, :, :) = Kkn .* ret(:, n, :, :);
end
end % function