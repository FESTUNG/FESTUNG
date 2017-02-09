function dataEval = evalFuncContAtEveryEdgeIntPoint(g, funcCont, Nlambda)
validateattributes(funcCont, {'function_handle'}, {}, mfilename, 'funcCont');

p = Nlambda-1;  qOrd = max(2*p, 1);  [Q1,~] = quadRule1D(qOrd);

% ord = max(ord,1);  [Q1, Q2, ~] = quadRule2D(ord);
% N = size(refElemPhiPhi, 1);

% K = g.numT;
% 
% dataEval = 
Kedge = g.numE;

dataEval = zeros( Kedge, size(Q1,2) );

for iE = 1:Kedge
    localIdx = g.E0E(iE,1);
    adjTri = g.T0E(iE, 1);
    
    [x1, x2] = gammaMap( localIdx, Q1 );
    
    F1 = @(X1, X2) g.B(adjTri,1,1)*X1 + g.B(adjTri,1,2)*X2 + g.coordV0T(adjTri,1,1)*ones(size(X1));
    F2 = @(X1, X2) g.B(adjTri,2,1)*X1 + g.B(adjTri,2,2)*X2 + g.coordV0T(adjTri,1,2)*ones(size(X1));
   
%     funcCont( F1(x1, x2), F2(x1, x2) )
    dataEval(iE, :)  = funcCont( F1(x1, x2), F2(x1, x2) ) ;
end
end % function