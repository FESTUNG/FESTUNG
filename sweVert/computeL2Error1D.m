function err = computeL2Error1D(g, dataDisc, funcCont, qOrd, basesOnQuad)
% Determine quadrature rule and physical coordinates
[Q, W] = quadRule1D(qOrd); X = g.mapRef2Phy(Q);
N = size(dataDisc,2);

% Evaluate analytical and discrete function
fExOnQuadPts = funcCont(X);
fApprxOnQuadPts = dataDisc * basesOnQuad.phi1D(:,1:N).'; % [K x R] = [K x N] * [N x R]

% Compute error
err = sqrt(dot((fApprxOnQuadPts - fExOnQuadPts).^2 * W.', g.areaT));
end  % function

