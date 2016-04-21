function ret = assembleMatEdgePhiIntNuPerQuad(g, markE0Tbdr, refEdgePhiIntPerQuad, areaNuE0Tbdr)
K = g.numT;
ret = cell(3,2);
for n = 1:3
	if nargin > 3
		ret{n,1} = kron(spdiags(areaNuE0Tbdr{n,1}, 0, K, K), refEdgePhiIntPerQuad(:,:,n));
		ret{n,2} = kron(spdiags(areaNuE0Tbdr{n,2}, 0, K, K), refEdgePhiIntPerQuad(:,:,n));
	else
		if isfield(g, 'areaNuE0T')
			ret{n,1} = kron(spdiags(g.areaNuE0T{n,1} .* markE0Tbdr(:,n), 0, K, K), refEdgePhiIntPerQuad(:,:,n));
			ret{n,2} = kron(spdiags(g.areaNuE0T{n,2} .* markE0Tbdr(:,n), 0, K, K), refEdgePhiIntPerQuad(:,:,n));
		else
			ret{n,1} = kron(spdiags(g.areaE0T(:,n) .* g.nuE0T(:,n,1) .* markE0Tbdr(:,n), 0, K, K), refEdgePhiIntPerQuad(:,:,n));
			ret{n,2} = kron(spdiags(g.areaE0T(:,n) .* g.nuE0T(:,n,2) .* markE0Tbdr(:,n), 0, K, K), refEdgePhiIntPerQuad(:,:,n));
		end % if
	end % if
end % for
end % function
