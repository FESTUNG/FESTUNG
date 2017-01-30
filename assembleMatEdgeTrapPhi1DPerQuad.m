function ret = assembleMatEdgeTrapPhi1DPerQuad(g, refEdgePhiIntPerQuad)
K = g.numT; barK = g.g1D.numT;
ret = cell(4,1);
for nn = 1 : 4
	ret{nn} = kron( spdiags(g.g1D.markT2DT.' * ( (g.markE0TE0T{nn} * ones(K,1)) .* g.areaE0T(:,nn) ), 0, barK, K), ...
                  refEdgePhiIntPerQuad(:,:,nn) );
end % for nn
end % function
