%
%
function ret = assembleMatEdgeMuPhiInt(g, markE0T, refEdgePhiIntMu)
K = g.numT;  N = size(refEdgePhiIntMu, 1);
Kedge = g.numE; Nlambda = size(refEdgePhiIntMu, 2);

% Check function arguments that are directly used
validateattributes(markE0T, {'logical'}, {'size', [K 3]}, mfilename, 'markE0T');
validateattributes(refEdgePhiIntMu, {'numeric'}, {'size', [N Nlambda 3 2]}, mfilename, 'refEdgePhiIntMu');

ret = sparse(K*N, Kedge*Nlambda);
for iE = 1:3
    Rkn1 = markE0T(:,iE) .*  g.flipArray(:, iE) .* g.areaE0T(:, iE);
    Rkn2 = markE0T(:,iE) .* ~g.flipArray(:, iE) .* g.areaE0T(:, iE);
    ret = ret + kron(sparse( 1:K, g.E0T(:, iE), Rkn1, K, Kedge ), ...
              refEdgePhiIntMu(:,:,iE, 1));
    ret = ret + kron(sparse( 1:K, g.E0T(:, iE), Rkn2, K, Kedge ), ...
              refEdgePhiIntMu(:,:,iE, 2));
end % for
end % function
