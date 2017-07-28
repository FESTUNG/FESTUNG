% TODO
function ret = assembleMatEdgePhiIntMuVal(g, markE0T, refEdgePhiIntMuPerQuad, valOnQuad)
% Extract dimensions
K = g.numT;
KEdge = g.numE;
[N, Nmu, ~, R] = size(refEdgePhiIntMuPerQuad{1});

% Check function arguments
validateattributes(markE0T, {'logical'}, {'size', [K 3]}, mfilename, 'markE0T');
validateattributes(refEdgePhiIntMuPerQuad, {'cell'}, {'size', [2 1]}, mfilename, 'refEdgePhiIntMuPerQuad');
validateattributes(refEdgePhiIntMuPerQuad{1}, {'numeric'}, {'size', [N Nmu 3 R]}, mfilename, 'refEdgePhiIntMuPerQuad{1}');
validateattributes(refEdgePhiIntMuPerQuad{2}, {'numeric'}, {'size', [N Nmu 3 R]}, mfilename, 'refEdgePhiIntMuPerQuad{2}');
validateattributes(valOnQuad, {'numeric'}, {'size', [K 3 R]}, mfilename, 'valOnQuad');

% Assemble matrix
ret = sparse(K*N, KEdge*Nmu);
for n = 1 : 3
  RknTimesVal = sparse(K*N, Nmu);
  markAreaE0T = markE0T(:, n) .* g.areaE0T(:, n);
  for l = 1 : 2
    markAreaSideE0T = markAreaE0T .*  g.markSideE0T(:, n, l);
    for r = 1 : R
      RknTimesVal = RknTimesVal + kron(markAreaSideE0T .* valOnQuad(:, n, r), refEdgePhiIntMuPerQuad{l}(:, :, n, r));
    end % for r
  end % for l
  ret = ret + kronVec(sparse(1:K, g.E0T(:, n), ones(K, 1), K, KEdge), RknTimesVal);
end % for n
end % function
