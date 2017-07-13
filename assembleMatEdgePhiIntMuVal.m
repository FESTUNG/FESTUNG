% TODO
function ret = assembleMatEdgePhiIntMuVal(g, markE0T, refEdgePhiIntMuOnQuad, valOnQuad)
% Extract dimensions
K = g.numT;
KEdge = g.numE;
[N, Nmu, ~, R] = size(refEdgePhiIntMuOnQuad{1});

% Check function arguments
validateattributes(markE0T, {'logical'}, {'size', [K 3]}, mfilename, 'markE0T');
validateattributes(refEdgePhiIntMuOnQuad, {'cell'}, {'size', [2 1]}, mfilename, 'refEdgePhiIntMuOnQuad');
validateattributes(refEdgePhiIntMuOnQuad{1}, {'numeric'}, {'size', [N Nmu 3 R]}, mfilename, 'refEdgePhiIntMuOnQuad{1}');
validateattributes(refEdgePhiIntMuOnQuad{2}, {'numeric'}, {'size', [N Nmu 3 R]}, mfilename, 'refEdgePhiIntMuOnQuad{2}');
validateattributes(valOnQuad, {'numeric'}, {'size', [K 3 R]}, mfilename, 'valOnQuad');

% Assemble matrix
ret = sparse(K*N, KEdge*Nmu);
for n = 1 : 3
  RknTimesVal = sparse(K*N, Nmu);
  markAreaE0T = markE0T(:, n) .* g.areaE0T(:, n);
  for l = 1 : 2
    markAreaSideE0T = markAreaE0T .*  g.markSideE0T(:, n, l);
    for r = 1 : R
      RknTimesVal = RknTimesVal + kron(markAreaSideE0T .* valOnQuad(:, n, r), refEdgePhiIntMuOnQuad{l}(:, :, n, r));
    end % for r
  end % for l
  ret = ret + kronVec(sparse(1:K, g.E0T(:, n), ones(K, 1), K, KEdge), RknTimesVal);
end % for n
end % function
