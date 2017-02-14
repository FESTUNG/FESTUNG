%
%
function ret = assembleMatEdgePhiIntMu(g, markE0T, refEdgePhiIntMu)
% Extract dimensions
K = g.numT;  N = size(refEdgePhiIntMu, 2);
Kedge = g.numE; Nlambda = size(refEdgePhiIntMu, 1);

% Check function arguments that are directly used
validateattributes(markE0T, {'logical'}, {'size', [K 3]}, mfilename, 'markE0T');
validateattributes(refEdgePhiIntMu, {'numeric'}, {'size', [Nlambda N 3 2]}, mfilename, 'refEdgePhiIntMu');

% Assemble matrix
ret = sparse(Kedge*Nlambda, K*N);
for iT = 1:K
    for iE = 1:3
        edgeNr = g.E0T(iT, iE);
        flip = 1;        
        if (g.T0E(edgeNr, 2) == iT)
           flip = 2; 
        end
        
        iTs = (iT-1)*N + 1;
        iTe = (iT)*N;
        iEs = (g.E0T(iT, iE) - 1)*Nlambda + 1;
        iEe = (g.E0T(iT, iE))*Nlambda;
        ret( iEs:iEe, iTs:iTe ) = ret( iEs:iEe, iTs:iTe ) + markE0T(iT,iE) .* g.areaE( edgeNr ) .* refEdgePhiIntMu(:,:, iE, flip);
    end
end


% for n = 1 : 3
%   ret = ret + kron(spdiags(markE0Tint(:,n),0,Kedge,K), refEdgePhiIntMu(:,:));
% end % for
end % function
