%
%
function ret = assembleMatEdgeMuPhiInt(g, refEdgePhiIntMu)
% Extract dimensions
K = g.numT;  N = size(refEdgePhiIntMu, 1);
KEdge = g.numE; Nlambda = size(refEdgePhiIntMu, 2);

% Check function arguments that are directly used
% validateattributes(markE0Tbdr, {'logical'}, {'size', [K 3]}, mfilename, 'markE0Tbdr');
validateattributes(refEdgePhiIntMu, {'numeric'}, {'size', [N Nlambda 3]}, mfilename, 'refEdgePhiIntMu');

% Assemble matrix
ret = sparse(K*N, KEdge*Nlambda);
% for n = 1 : 3
%   %ret = ret + stab .* ( kron(spdiags(markE0Tbdr(:,n),0,K,K), refEdgePhiIntMu(:,:)) );
%   ret = ret + ( kron(spdiags(g.areaE, 0, K, KEdge), refEdgePhiIntMu(:,:, n) ) );
%   
% end % for
for iT = 1:K
    for iE = 1:3
        iTs = (iT-1)*N + 1;
        iTe = (iT)*N;
        iEs = (g.E0T(iT, iE) - 1)*Nlambda + 1;
        iEe = (g.E0T(iT, iE))*Nlambda;
        ret( iTs:iTe,  iEs:iEe  ) = ret( iTs:iTe, iEs:iEe ) +  refEdgePhiIntMu(:,:, iE);
    end
end

end % function
