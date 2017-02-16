function dataEval = evalFluxContAtEveryEdgeIntPoint(g, problemData, fluxCont, cCont, Nlambda)
validateattributes(fluxCont, {'function_handle'}, {}, mfilename, 'fluxCont');

p = Nlambda-1;  qOrd = max(2*p, 1);  [Q,~] = quadRule1D(qOrd);

% ord = max(ord,1);  [Q1, Q2, ~] = quadRule2D(ord);
% N = size(refElemPhiPhi, 1);

K = g.numT;
%
% dataEval =
Kedge = g.numE;

% dataEval = zeros( 2, size(Q,2), K, 3 );

dataEval = zeros( K, 3, 2, size(Q,2) );


for iE = 1:3
    for iT = 1:K
        edgeNr = g.E0T(iT, iE);
        
        localIdx = iE;
        adjTri = iT;
        
        %         localIdx = g.E0E(iE,1);
        %         adjTri = g.T0E(iE, 1);
        
        if ( g.markE0TbdrD(iT, iE) )
            
            [x1, x2] = gammaMap( localIdx, Q );
            
            F1 = @(X1, X2) g.B(adjTri,1,1)*X1 + g.B(adjTri,1,2)*X2 + g.coordV0T(adjTri,1,1)*ones(size(X1));
            F2 = @(X1, X2) g.B(adjTri,2,1)*X1 + g.B(adjTri,2,2)*X2 + g.coordV0T(adjTri,1,2)*ones(size(X1));
            
            X1D = F1(x1, x2) ;
            X2D = F2(x1, x2) ;
            test = fluxCont( F1(x1, x2), F2(x1, x2), 1 );
            %     cCont(iE, :)
            %     fluxCont( F1(x1, x2), F1(x1, x2), cCont(iE, :) )
%             dataEval(:, :, iT, iE)  = ~g.markE0Tint(iT, iE) * fluxCont( F1(x1, x2), F2(x1, x2), cCont(iE, :) );
            dataEval( iT, iE, :, :) = ~g.markE0Tint(iT, iE) * fluxCont( F1(x1, x2), F2(x1, x2), cCont(iE, :) );
        else
            %         if ( g.markE0Tint(iT, iE) )
            [x1, x2] = gammaMap( localIdx, Q );
            
            F1 = @(X1, X2) g.B(adjTri,1,1)*X1 + g.B(adjTri,1,2)*X2 + g.coordV0T(adjTri,1,1)*ones(size(X1));
            F2 = @(X1, X2) g.B(adjTri,2,1)*X1 + g.B(adjTri,2,2)*X2 + g.coordV0T(adjTri,1,2)*ones(size(X1));
            
            bases = problemData.basesOnGamma.phi1D{qOrd};
            cDiscLambda = problemData.cDiscLambda;
            
            lambdaLocal = cDiscLambda(edgeNr,:) * bases';
            
%             dataEval(:, :, iT, iE) = fluxCont( F1(x1, x2), F2(x1, x2), lambdaLocal );
            dataEval( iT, iE, :, :) = fluxCont( F1(x1, x2), F2(x1, x2), lambdaLocal );

        end
        
        %     dataEval(iE, :, : )  = fluxCont( F1(x1, x2), F2(x1, x2), cCont(iE, :) );
    end
end

% dataEvalTest = zeros( 2, size(Q,2), K, 3 );
% F1 = @(X1, X2) g.B(:,1,1)*X1 + g.B(:,1,2)*X2 + g.coordV0T(:,1,1)*ones(size(X1));
% F2 = @(X1, X2) g.B(:,2,1)*X1 + g.B(:,2,2)*X2 + g.coordV0T(:,1,2)*ones(size(X1));
% for n = 1:3
%     [Q1, Q2] = gammaMap(n, Q);
%     Kkn = ~g.markE0Tint(:, n) .* g.markE0TbdrD(:, n);
% %     dataEvalTest(:, :, :, n) = reshape( fluxCont( F1(Q1, Q2), F2(Q1, Q2), cCont(iE, :) ), 2, 2, K);
%   dataEvalTest(:, :, :, n) = reshape( kron( ones(2,2), Kkn ) .* fluxCont( F1(Q1, Q2), F2(Q1, Q2), cCont(:, :) ), 2, 2, K);
% 
% end
% 
% lambdaEval = cDiscLambda(:,:) * bases';
% for n = 1:3
%     [Q1, Q2] = gammaMap(n, Q);
%     Kkn = ~g.markE0TbdrD(:, n);
% %     dataEvalTest(:, :, :, n) = Kkn * fluxCont( F1(Q1, Q2), F2(Q1, Q2), lambdaEval );
%     dataEvalTest(:, :, :, n) = reshape( kron( ones(2,2), Kkn ) .* fluxCont( F1(Q1, Q2), F2(Q1, Q2), lambdaEval(n,:) ), 2, 2, K);
% end


end % function