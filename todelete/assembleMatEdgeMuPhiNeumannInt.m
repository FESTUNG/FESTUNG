%
%
function ret = assembleMatEdgeMuPhiNeumannInt(g, markE0TbdrN, refEdgePhiIntMu)
% Extract dimensions
K = g.numT;  N = size(refEdgePhiIntMu, 1);
KEdge = g.numE; Nlambda = size(refEdgePhiIntMu, 2);

% Check function arguments that are directly used
% validateattributes(markE0Tbdr, {'logical'}, {'size', [K 3]}, mfilename, 'markE0Tbdr');
validateattributes(refEdgePhiIntMu, {'numeric'}, {'size', [N Nlambda 3 2]}, mfilename, 'refEdgePhiIntMu');

% warning('area removed from assembleMatEdgeMuPhiInt');

% Assemble matrix
ret = sparse(K*N, KEdge*Nlambda);
for iT = 1:K
    for iE = 1:3
        edgeNr = g.E0T(iT, iE);
        flip = 1;
        if (g.T0E(edgeNr, 2) == iT)
            flip = 2;
        end
        iTs = (iT-1)*N + 1;
        iTe = (iT)*N;
        iEs = (edgeNr - 1)*Nlambda + 1;
        iEe = (edgeNr)*Nlambda;

        ret( iTs:iTe,  iEs:iEe  ) = ret( iTs:iTe, iEs:iEe ) + g.areaE0T( iT, iE ) .* markE0TbdrN(iT, iE) .* refEdgePhiIntMu(:,:, iE, flip);
    end
end


end % function
