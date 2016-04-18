function ret = assembleMatEdgePhiIntFuncDiscConstNu(g, markE0Tbdr, refEdgePhiIntPerQuad, areaNuE0Tbdr)
K = g.numT;
N = size(refEdgePhiIntPerQuad, 1);
ret = cell(2,1);
ret{1} = sparse(K*N,1); ret{2} = sparse(K*N,1);
for n = 1:3
	if nargin > 3
		ret{1} = ret{1} + kron(areaNuE0Tbdr{n,1}, refEdgePhiIntPerQuad(:,:,n));
		ret{2} = ret{2} + kron(areaNuE0Tbdr{n,2}, refEdgePhiIntPerQuad(:,:,n));
	else
		ret{1} = ret{1} + kron(g.areaE0T(:,n) .* g.nuE0T(:,n,1) .* markE0Tbdr(:,n), refEdgePhiIntPerQuad(:,:,n));
		ret{2} = ret{2} + kron(g.areaE0T(:,n) .* g.nuE0T(:,n,2) .* markE0Tbdr(:,n), refEdgePhiIntPerQuad(:,:,n));
	end % if
end % for
end % function
