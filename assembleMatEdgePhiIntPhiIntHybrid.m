%
function ret = assembleMatEdgePhiIntPhiIntHybrid(g, refEdgePhiIntPhiInt)
% Extract dimensions
K = g.numT;  N = size(refEdgePhiIntPhiInt, 1);

% Check function arguments that are directly used
validateattributes(refEdgePhiIntPhiInt, {'numeric'}, {'size', [N N 3]}, mfilename, 'refEdgePhiIntPhiInt');

ret = sparse(K*N, K*N);
for n = 1 : 3
    ret = ret + kron(spdiags( g.markE0Tint(:, n) .* g.areaE( g.E0T(:, n) ) ,0,K,K), refEdgePhiIntPhiInt(:,:,n) );
end % for

end % function
