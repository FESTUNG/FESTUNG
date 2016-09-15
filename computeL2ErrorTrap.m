function err = computeL2ErrorTrap(g, dataDisc, funcCont, qOrd, basesOnQuad)
% Determine quadrature rule and physical coordinates
[Q, W] = quadRule1D(qOrd); [Q1, Q2] = meshgrid(Q); W = W' * W;
Q1 = Q1(:)'; Q2 = Q2(:)'; W = W(:)';
X1 = g.mapRef2Phy(1,Q1,Q2); X2 = g.mapRef2Phy(2,Q1,Q2);
N = size(dataDisc,2);

% Evaluate analytical and discrete function
fExOnQuadPts = funcCont(X1, X2);
fApprxOnQuadPts = dataDisc * basesOnQuad.phi2D(:,1:N).'; % [K x R] = [K x N] * [N x R]

% Compute error
err = sqrt(dot((fApprxOnQuadPts - fExOnQuadPts).^2 * W.', g.areaT));
end  % function

