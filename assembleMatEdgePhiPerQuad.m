function ret = assembleMatEdgePhiPerQuad(g, refEdgePhiIntPerQuad)
K = g.numT;
ret = cell(3,3);
for nn = 1 : 3
	for np = 1: 3
		if isfield(g, 'areaMarkE0T')
			ret{nn,np} = 0.5 * kron( spdiags(g.areaMarkE0T{nn,np}, 0, K, K), refEdgePhiIntPerQuad(:,:,nn) );
		else
			ret{nn,np} = 0.5 * kron( spdiags((g.markE0TE0T{nn,np} * ones(K,1)) .* g.areaE0T(:,nn), 0, K, K), refEdgePhiIntPerQuad(:,:,nn) );
		end % if
	end % for
end % for
end % function
