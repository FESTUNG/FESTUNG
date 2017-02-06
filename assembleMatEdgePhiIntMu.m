%
%
function ret = assembleMatEdgePhiIntMu(g, markE0Tbdr, refEdgePhiIntMu)
% Extract dimensions
K = g.numT;  N = size(refEdgePhiIntMu, 2);

% Check function arguments that are directly used
validateattributes(markE0Tbdr, {'logical'}, {'size', [K 3]}, mfilename, 'markE0Tbdr');
validateattributes(refEdgePhiIntMu, {'numeric'}, {'size', [N N 3]}, mfilename, 'refEdgePhiIntMu');

% Assemble matrix
ret = sparse(K*N, K*N);
for n = 1 : 3
  ret = ret + kron(spdiags(markE0Tbdr(:,n),0,K,K), refEdgePhiIntMu(:,:));
end % for
end % function
