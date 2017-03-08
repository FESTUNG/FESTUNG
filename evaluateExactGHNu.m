function rhs = evaluateExactGHNu(g, t, gravity, markE0T, funcCont, qOrd, basesOnQuad)

validateattributes(funcCont, {'function_handle'}, {}, mfilename, 'funcCont');
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad');
[Q, W] = quadRule1D(max(qOrd,1));
rhs = zeros(size(markE0T, 1), size(basesOnQuad.phi2D, 2));
for i = [3, 4]
    switch i
        case 3
            Q1 = g.mapRef2Phy(1,ones(size(Q)),Q);
        case 4
            Q1 = g.mapRef2Phy(1,zeros(size(Q)),Q);
    end
    rhs = rhs + gravity * funcCont(t,Q1) * (repmat(W.', 1, size(basesOnQuad.phi1D(:,:,i),2)) ...
            .* basesOnQuad.phi1D(:,:,i)) .* kron((markE0T(:,i) .* g.areaE0T(:,i) .* g.nuE0T(:,i,1)), ones(1, size(basesOnQuad.phi2D, 2)));
end

rhs = reshape(rhs.', size(rhs,1) * size(rhs,2), 1);

end