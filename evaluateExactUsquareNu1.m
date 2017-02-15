function [ rhs ] = evaluateExactUsquareNu1(g, t, markE0T, funcCont, qOrd, basesOnQuad)

validateattributes(funcCont, {'function_handle'}, {}, mfilename, 'funcCont');
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad');
[Q, W] = quadRule1D(max(qOrd,1));
rhs = zeros(size(markE0T, 1), size(basesOnQuad.phi2D, 2));
for i = 1 : 4
    switch i
        case 1
            Q1 = g.mapRef2Phy(1,Q,zeros(size(Q)));
            Q2 = g.mapRef2Phy(2,Q,zeros(size(Q)));
        case 2
            Q1 = g.mapRef2Phy(1,Q,ones(size(Q)));
            Q2 = g.mapRef2Phy(2,Q,ones(size(Q)));
        case 3
            Q1 = g.mapRef2Phy(1,ones(size(Q)),Q);
            Q2 = g.mapRef2Phy(2,ones(size(Q)),Q);
        case 4
            Q1 = g.mapRef2Phy(1,zeros(size(Q)),Q);
            Q2 = g.mapRef2Phy(2,zeros(size(Q)),Q);
    end
    rhs = rhs + funcCont(t,Q1,Q2).^2 * (repmat(W.', 1, size(basesOnQuad.phi1D(:,:,i),2)) ...
            .* basesOnQuad.phi1D(:,:,i)) .* kron((markE0T(:,i) .* g.areaE0T(:,i) .* g.nuE0T(:,i,1)), ones(1, size(basesOnQuad.phi2D, 2)));
end

rhs = reshape(rhs.', size(rhs,1) * size(rhs,2), 1);

end