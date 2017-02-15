function dataEval = evalFluxContAtEveryEdgeIntPoint(g, problemData, fluxCont, cCont, Nlambda)
validateattributes(fluxCont, {'function_handle'}, {}, mfilename, 'fluxCont');

p = Nlambda-1;  qOrd = max(2*p, 1);  [Q1,~] = quadRule1D(qOrd);

% ord = max(ord,1);  [Q1, Q2, ~] = quadRule2D(ord);
% N = size(refElemPhiPhi, 1);

K = g.numT;
%
% dataEval =
Kedge = g.numE;

dataEval = zeros( 2, size(Q1,2), K, 3 );

for iT = 1:K
    for iE = 1:3
        edgeNr = g.E0T(iT, iE);
        
        localIdx = iE;
        adjTri = iT;
        
        %         localIdx = g.E0E(iE,1);
        %         adjTri = g.T0E(iE, 1);
        
        if ( g.markE0TbdrD(iT, iE) )
            
            [x1, x2] = gammaMap( localIdx, Q1 );
            
            F1 = @(X1, X2) g.B(adjTri,1,1)*X1 + g.B(adjTri,1,2)*X2 + g.coordV0T(adjTri,1,1)*ones(size(X1));
            F2 = @(X1, X2) g.B(adjTri,2,1)*X1 + g.B(adjTri,2,2)*X2 + g.coordV0T(adjTri,1,2)*ones(size(X1));
            
            X1D = F1(x1, x2) ;
            X2D = F2(x1, x2) ;
            test = fluxCont( F1(x1, x2), F2(x1, x2), 1 );
            %     cCont(iE, :)
            %     fluxCont( F1(x1, x2), F1(x1, x2), cCont(iE, :) )
            dataEval(:, :, iT, iE)  = ~g.markE0Tint(iT, iE) * fluxCont( F1(x1, x2), F2(x1, x2), cCont(iE, :) );
        else
            %         if ( g.markE0Tint(iT, iE) )
            [x1, x2] = gammaMap( localIdx, Q1 );
            
            F1 = @(X1, X2) g.B(adjTri,1,1)*X1 + g.B(adjTri,1,2)*X2 + g.coordV0T(adjTri,1,1)*ones(size(X1));
            F2 = @(X1, X2) g.B(adjTri,2,1)*X1 + g.B(adjTri,2,2)*X2 + g.coordV0T(adjTri,1,2)*ones(size(X1));
            
            bases = problemData.basesOnGamma.phi1D{qOrd};
            cDiscLambda = problemData.cDiscLambda;
            
            lambdaLocal = cDiscLambda(edgeNr,:) * bases';
            
            dataEval(:, :, iT, iE) = fluxCont( F1(x1, x2), F2(x1, x2), lambdaLocal );
        end
        
        %     dataEval(iE, :, : )  = fluxCont( F1(x1, x2), F2(x1, x2), cCont(iE, :) );
    end
end
end % function