% TODO Check this
% Compute global mas matrix of lambda
% If it is an internal edge, the entry has to be scaled by 2*alpha
% check the element transformation
function ret = assembleMatEdgeMuMu(g, alpha, refEdgeMuMu)
KEdge = g.numE;
ret = kron(spdiags(g.areaE, 0, KEdge, KEdge), refEdgeMuMu);
end % function
