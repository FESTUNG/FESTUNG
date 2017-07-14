% TODO
function dataDisc = projectFuncCont2EdgeDataDisc(g, funcCont, ord, refEdgePhiPhi, basesOnGamma)
validateattributes(funcCont, {'function_handle'}, {}, mfilename, 'funcCont');
validateattributes(basesOnGamma, {'struct'}, {}, mfilename, 'basesOnGamma');
ord = max(ord,1);  [Q, W] = quadRule1D(ord);
N = size(refEdgePhiPhi, 1);
rhs = zeros( g.numE, N );
for n = 1:3
    [Q1, Q2] = gammaMap(n, Q);
    rhs(g.E0T(:, n),:) = repmat( W, g.numT, 1) .* funcCont(g.mapRef2Phy(1, Q1, Q2), g.mapRef2Phy(2, Q1, Q2) ) * basesOnGamma.phi1D{ord};    
end
dataDisc = rhs / refEdgePhiPhi;
end % function
