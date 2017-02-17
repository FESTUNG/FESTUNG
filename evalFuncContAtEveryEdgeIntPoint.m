function dataEval = evalFuncContAtEveryEdgeIntPoint(g, funcCont, Nlambda)
validateattributes(funcCont, {'function_handle'}, {}, mfilename, 'funcCont');

p = Nlambda-1;  qOrd = max(2*p, 1);  [Q1,~] = quadRule1D(qOrd);

% ord = max(ord,1);  [Q1, Q2, ~] = quadRule2D(ord);
% N = size(refElemPhiPhi, 1);

% K = g.numT;
% 
% dataEval = 
K = g.numT;
Kedge = g.numE;

% dataEval = zeros( Kedge, size(Q1,2) );
% 
% for iE = 1:Kedge
%     localIdx = g.E0E(iE,1);
%     adjTri = g.T0E(iE, 1);
%     
%     [x1, x2] = gammaMap( localIdx, Q1 );
%     
%     F1 = @(X1, X2) g.B(adjTri,1,1)*X1 + g.B(adjTri,1,2)*X2 + g.coordV0T(adjTri,1,1)*ones(size(X1));
%     F2 = @(X1, X2) g.B(adjTri,2,1)*X1 + g.B(adjTri,2,2)*X2 + g.coordV0T(adjTri,1,2)*ones(size(X1));
%    
% %     funcCont( F1(x1, x2), F2(x1, x2) )
%     dataEval(iE, :)  = funcCont( F1(x1, x2), F2(x1, x2) ) ;
% end

dataEval = zeros( K, size(Q1,2), 3 );

F1 = @(X1, X2) g.B(:,1,1)*X1 + g.B(:,1,2)*X2 + g.coordV0T(:,1,1)*ones(size(X1));
F2 = @(X1, X2) g.B(:,2,1)*X1 + g.B(:,2,2)*X2 + g.coordV0T(:,1,2)*ones(size(X1));
for n = 1:3
    [x1, x2] = gammaMap( n, Q1 );
   
%     funcCont( F1(x1, x2), F2(x1, x2) )
    dataEval(:, :, n)  = funcCont( F1(x1, x2), F2(x1, x2) ) ;
end


% dataEvalTest = zeros( Kedge, size(Q1,2) );
% for n = 1:3
%     [x1, x2] = gammaMap( n, Q1 );
%     F1 = @(X1, X2) g.B(:,1,1)*X1 + g.B(:,1,2)*X2 + g.coordV0T(:,1,1)*ones(size(X1));
%     F2 = @(X1, X2) g.B(:,2,1)*X1 + g.B(:,2,2)*X2 + g.coordV0T(:,1,2)*ones(size(X1));
%     
% %     dataEvalTest = funcCont( F1(x1, x2), F2(x1, x2) )
%     
%     dataEvalTest(:, :) = dataEvalTest(:, :) + sparse( g.E0T(:, n), 1:K, ones(K,1), Kedge, K) * funcCont( F1(x1, x2), F2(x1, x2) );
% end

% dataEvalTest = zeros( K, 3, size(Q1,2) );
% for n = 1:3
%     [x1, x2] = gammaMap( n, Q1 );
%     F1 = @(X1, X2) g.B(:,1,1)*X1 + g.B(:,1,2)*X2 + g.coordV0T(:,1,1)*ones(size(X1));
%     F2 = @(X1, X2) g.B(:,2,1)*X1 + g.B(:,2,2)*X2 + g.coordV0T(:,1,2)*ones(size(X1));
%     
% %     dataEvalTest = funcCont( F1(x1, x2), F2(x1, x2) )
%     
% %     dataEvalTest(:, :) = dataEvalTest(:, :) + sparse( g.E0T(:, n), 1:K, ones(K,1), Kedge, K) * funcCont( F1(x1, x2), F2(x1, x2) );
%     dataEvalTest(:, n, :) =  funcCont( F1(x1, x2), F2(x1, x2) );
% end

end % function