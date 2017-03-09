function ret = assembleMatEdgeTrapPhiPerQuad(g, refEdgePhiIntPerQuad)
K = g.numT;
ret = cell(4,1);
for n = 1 : 4
	ret{n} = kron( spdiags((g.markE0TE0T{n} * ones(K,1)) .* g.areaE0T(:,n), 0, K, K), refEdgePhiIntPerQuad(:,:,n) );
end % for nn
end % function
