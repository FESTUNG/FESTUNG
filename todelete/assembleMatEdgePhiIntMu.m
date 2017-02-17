%
%
function ret = assembleMatEdgePhiIntMu(g, markE0T, refEdgePhiIntMu)
% Extract dimensions
K = g.numT;  N = size(refEdgePhiIntMu, 2);
Kedge = g.numE; Nlambda = size(refEdgePhiIntMu, 1);

% Check function arguments that are directly used
validateattributes(markE0T, {'logical'}, {'size', [K 3]}, mfilename, 'markE0T');
validateattributes(refEdgePhiIntMu, {'numeric'}, {'size', [Nlambda N 3 2]}, mfilename, 'refEdgePhiIntMu');

ret = sparse(Kedge*Nlambda, K*N);
for iE = 1:3
    Rkn1 = markE0T(:,iE) .*  g.flipArray(:, iE) .* g.areaE0T(:, iE);
    Rkn2 = markE0T(:,iE) .* ~g.flipArray(:, iE) .* g.areaE0T(:, iE);
    ret = ret + kron(sparse( g.E0T(:, iE), 1:K, Rkn1, Kedge, K ), ...
              refEdgePhiIntMu(:,:,iE, 1));
    ret = ret + kron(sparse( g.E0T(:, iE), 1:K, Rkn2, Kedge, K ), ...
              refEdgePhiIntMu(:,:,iE, 2));
end % for
end % function
