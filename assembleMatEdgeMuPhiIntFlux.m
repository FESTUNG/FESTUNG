function ret = assembleMatEdgeMuPhiIntFlux( g, markE0Tbdr, N, Nlambda, uEdge, Sbar )

%Assert that the length of uEval and hatSbarOnQuad is ok
K = g.numT;
KEdge = g.numE;

% Assemble matrix
ret = cell(2,1);
ret{1} = sparse( K*N, KEdge*Nlambda);
ret{2} = sparse( K*N, KEdge*Nlambda);

for iT = 1:K
    for iE = 1:3
        edgeNr = g.E0T(iT, iE);
        
        uEdgeLocal = uEdge( :, :, iT, iE);
        flip = 1;
        %         if (g.T0E(edgeNr, 2) == iT)
        if (g.T0E(edgeNr, 2) == iT)
            %         if (2 == iT)
            flip = 2;
%             uEdgeLocal = fliplr( uEdgeLocal );
        end
        
        iTs = (iT-1)*N + 1;
        iTe = (iT)*N;
        iEs = (g.E0T(iT, iE) - 1)*Nlambda + 1;
        iEe = (g.E0T(iT, iE))*Nlambda;
        
        for iDim = 1:2
            tmp = zeros(N, Nlambda);
            for i = 1:N
                for j=1:Nlambda
                    %                     tmp(i,j) = uEdge( iDim, :, iT, iE) * Sbar( :, i, j, iDim, flip);
                    tmp(i,j) = uEdgeLocal( iDim, :) .* g.nuE0T( iT, iE , iDim) * Sbar( :, i, j, iE, flip);
                end
            end
            
            ret{iDim}( iTs:iTe,  iEs:iEe ) = ret{iDim}( iTs:iTe, iEs:iEe ) ...
                + markE0Tbdr( iT, iE ) .* g.areaE( edgeNr ) .* tmp;
        end
    end
end
end