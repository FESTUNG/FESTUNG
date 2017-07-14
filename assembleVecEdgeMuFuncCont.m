% TODO
function ret = assembleVecEdgeMuFuncCont(g, markE0T, funcCont, basesOnQuad, qOrd)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
validateattributes(funcCont, {'function_handle'}, {}, mfilename, 'funcCont');
validateattributes(markE0T, {'logical'}, {'size', [g.numT 3]}, mfilename, 'markE0T');

% Determine quadrature rule
[Q, W] = quadRule1D(qOrd);
[R, Nmu] = size(basesOnQuad.mu{qOrd});

% Assemble vector
ret = zeros(g.numE, Nmu);
for n = 1 : 3
  [Q1, Q2] = gammaMap(n, Q);
  funcOnQuad = funcCont(g.mapRef2Phy(1, Q1, Q2), g.mapRef2Phy(2, Q1, Q2));
  Kkn = markE0T(:, n) .* g.areaE0T(:,n);
  ret(g.E0T(:, n), :) = ret(g.E0T(:, n), :) + (repmat(Kkn, 1, R) .* funcOnQuad) ...
                          * ( repmat(W', 1, Nmu) .* basesOnQuad.mu{qOrd} );
end % for

ret = reshape(ret', [], 1);
end


