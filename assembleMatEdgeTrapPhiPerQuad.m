function ret = assembleMatEdgeTrapPhiPerQuad(g, refEdgePhiIntPerQuad)
K = g.numT;
ret = cell(4,1);
for nn = 1 : 4
	ret{nn} = kron( spdiags((g.markE0TE0T{nn} * ones(K,1)) .* g.areaE0T(:,nn), 0, K, K), refEdgePhiIntPerQuad(:,:,nn) );
end % for nn
end % function
