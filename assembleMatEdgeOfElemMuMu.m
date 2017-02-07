% TODO Check this
% Compute global mas matrix of lambda
% If it is an internal edge, the entry has to be scaled by 2*alpha
% check the element transformation
function ret = assembleMatEdgeOfElemMuMu(g, alpha, refEdgeMuMu)

K = g.numT;  Nlambda = size(refEdgeMuMu, 1);
KEdge = g.numE;
ret = sparse(KEdge*Nlambda, KEdge*Nlambda);

%Interior edges
for iT = 1:K
    for iE = 1:3
%         iTs = (iT-1)*N + 1;
%         iTe = (iT)*N;
        iEs = (g.E0T(iT, iE) - 1)*Nlambda + 1;
        iEe = (g.E0T(iT, iE))*Nlambda;
        ret( iEs:iEe,  iEs:iEe  ) = ret( iEs:iEe,  iEs:iEe ) + g.areaE( g.E0T(iT, iE) ) .* ( g.markE0Tint(iT, iE) ) .* alpha .* refEdgeMuMu;        
    end
end

%Exterior edges
for iT = 1:K
    for iE = 1:3
%         iTs = (iT-1)*N + 1;
%         iTe = (iT)*N;
        iEs = (g.E0T(iT, iE) - 1)*Nlambda + 1;
        iEe = (g.E0T(iT, iE))*Nlambda;
        ret( iEs:iEe,  iEs:iEe  ) = ret( iEs:iEe,  iEs:iEe ) + g.areaE( g.E0T(iT, iE) ) .* ( ~g.markE0Tint(iT, iE) ) .* refEdgeMuMu;        
    end
end

% ret = kron(spdiags(g.areaE, 0, KEdge, KEdge), refEdgeMuMu);
end % function
