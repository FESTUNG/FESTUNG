% TODO
function dataDisc = projectFuncCont2EdgeDataDisc(g, funcCont, qOrd, refEdgePhiPhi, basesOnQuad)
validateattributes(funcCont, {'function_handle'}, {}, mfilename, 'funcCont');
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnGamma');
qOrd = max(qOrd,1);  [Q, W] = quadRule1D(qOrd);
N = size(refEdgePhiPhi, 1);
rhs = zeros(g.numE, N);
for n = 1:3
  [Q1, Q2] = gammaMap(n, Q);
  rhs(g.E0T(:, n),:) = repmat(W, g.numT, 1) .* funcCont(g.mapRef2Phy(1, Q1, Q2), g.mapRef2Phy(2, Q1, Q2) ) * basesOnQuad.mu{qOrd};
end
dataDisc = rhs / refEdgePhiPhi;
end % function
