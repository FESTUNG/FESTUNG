function ret = assembleMatEdgeMuPhiIntFlux( g, markE0Tbdr,uEdge, Sbar )

%Assert that the length of uEval and hatSbarOnQuad is ok
K = g.numT;
KEdge = g.numE;

[N, Nlambda, ~, R] = size(Sbar{1});
% Assemble matrix
ret = cell(2,1);
ret{1} = sparse( K*N, KEdge*Nlambda);
ret{2} = sparse( K*N, KEdge*Nlambda);
for iDim = 1:2
    for iE = 1:3
        Rkn1 = markE0Tbdr(:,iE) .* g.areaE0T(:, iE) .*  g.flipArray(:, iE) .* g.nuE0T( :, iE , iDim);
        Rkn2 = markE0Tbdr(:,iE) .* g.areaE0T(:, iE) .* ~g.flipArray(:, iE) .* g.nuE0T( :, iE , iDim);
        for r = 1:R
            ret{iDim} = ret{iDim} + ...
                kron( sparse( 1:K, g.E0T(:, iE), ones(K, 1) .* uEdge{iDim}( r, :, iE)' .* Rkn1 , K, KEdge ), ...
                Sbar{1}( :, :, iE, r) );
            ret{iDim} = ret{iDim} + ...
                kron( sparse( 1:K, g.E0T(:, iE), ones(K, 1) .* uEdge{iDim}( r, :, iE)' .* Rkn2 , K, KEdge ), ...
                Sbar{2}( :, :, iE, r) );
        end
    end
end
end