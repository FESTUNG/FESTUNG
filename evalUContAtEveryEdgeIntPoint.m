function dataEval = evalUContAtEveryEdgeIntPoint(g, cCont, Nlambda)
validateattributes(cCont, {'function_handle'}, {}, mfilename, 'cCont');

p = Nlambda-1;  qOrd = max(2*p, 1);  [Q1,~] = quadRule1D(qOrd);
K = g.numT;
Kedge = g.numE;

% dataEval = zeros( 2, size(Q1,2), K, 3 );
dataEval = cell(2,1);
dataEval{1} = zeros( size(Q1,2), K, 3 );
dataEval{2} = zeros( size(Q1,2), K, 3 );


for iT = 1:K
    for iE = 1:3
        localIdx = iE;
        adjTri = iT;

        [x1, x2] = gammaMap( localIdx, Q1 );
        
        F1 = @(X1, X2) g.B(adjTri,1,1)*X1 + g.B(adjTri,1,2)*X2 + g.coordV0T(adjTri,1,1)*ones(size(X1));
        F2 = @(X1, X2) g.B(adjTri,2,1)*X1 + g.B(adjTri,2,2)*X2 + g.coordV0T(adjTri,1,2)*ones(size(X1));
        
%         dataEval(:, :, iT, iE)  = cCont( F1(x1, x2), F2(x1, x2) );
        res = cCont( F1(x1, x2), F2(x1, x2) );
        dataEval{1}(:, iT, iE) = res(1, :);
        dataEval{2}(:, iT, iE) = res(2, :);

%         dataEval(:, :, iT, iE)  = (~g.markE0Tint(iT, iE)) * cCont( F1(x1, x2), F2(x1, x2) );
    end
end

% for iT = 1:K
%     for iE = 1:3
%         localIdx = iE;
%         adjTri = iT;
% 
% %         [x1, x2] = gammaMap( localIdx, Q1 );
% %         
% %         F1 = @(X1, X2) g.B(adjTri,1,1)*X1 + g.B(adjTri,1,2)*X2 + g.coordV0T(adjTri,1,1)*ones(size(X1));
% %         F2 = @(X1, X2) g.B(adjTri,2,1)*X1 + g.B(adjTri,2,2)*X2 + g.coordV0T(adjTri,1,2)*ones(size(X1));
%         
%         
% 
%         dataEval(:, :, iT, iE) = (~g.markE0Tint(iT, iE)) * cCont( F1(x1, x2), F2(x1, x2) );
%     end
% end

% for iE = 1:Kedge
%     localIdx = g.E0E(iE,1);
%     adjTri = g.T0E(iE, 1);
%
%     [x1, x2] = gammaMap( localIdx, Q1 );
%
%     F1 = @(X1, X2) g.B(adjTri,1,1)*X1 + g.B(adjTri,1,2)*X2 + g.coordV0T(adjTri,1,1)*ones(size(X1));
%     F2 = @(X1, X2) g.B(adjTri,2,1)*X1 + g.B(adjTri,2,2)*X2 + g.coordV0T(adjTri,1,2)*ones(size(X1));
%
%     dataEval(:, :, iE)  = cCont( F1(x1, x2), F2(x1, x2) );
% end
end % function