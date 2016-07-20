function cDG = projectAlg2DG(g, cAlg, ord, hatM, basesOnQuad)
% global gPhi2D
ord = max(ord,1);  [Q1, Q2, W] = quadRule2D(ord);
N = size(hatM, 1);
F1 = @(X1, X2) g.B(:,1,1)*X1 + g.B(:,1,2)*X2 + g.coordV0T(:,1,1)*ones(size(X1));
F2 = @(X1, X2) g.B(:,2,1)*X1 + g.B(:,2,2)*X2 + g.coordV0T(:,1,2)*ones(size(X1));
rhs = cAlg(F1(Q1, Q2), F2(Q1, Q2)) * (repmat(W.', 1, N) .* basesOnQuad.phi2D{ord});
cDG = rhs / hatM;
end % function
