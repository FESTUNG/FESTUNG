% TODO Check this
%
function ret = assembleMatEdgePhiPhiHybrid(g, refEdgePhiPhi)
KEdge = g.numE;
ret = kron(spdiags(g.areaE, 0, KEdge, KEdge), refEdgePhiPhi);
end % function
