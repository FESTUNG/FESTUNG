function ret = assembleMatEdgeMuPhiIntFlux( g, N, Nlambda, uEval, Sbar )

%Assert that the length of uEval and hatSbarOnQuad is ok
K = g.numT;
KEdge = g.numE;

% Assemble matrix
ret = cell(2,1); 
ret{1} = sparse( K*N, KEdge*Nlambda);
ret{2} = sparse( K*N, KEdge*Nlambda);

for iT = 1:K
    for iE = 1:3
        iTs = (iT-1)*N + 1;
        iTe = (iT)*N;
        iEs = (g.E0T(iT, iE) - 1)*Nlambda + 1;
        iEe = (g.E0T(iT, iE))*Nlambda;
        
        for iDim = 1:2
            
            tmp = zeros(N, Nlambda);
            for i = 1:N
                for j=1:Nlambda
                    tmp(i,j) = uEval( iT, :, iDim)  * Sbar( :, i, j, iDim) ;
                end
            end
            
            
            ret{iDim}(iTs:iTe,  iEs:iEe  ) = ret{iDim}( iTs:iTe, iEs:iEe ) +  tmp;
        end
    end
end


end