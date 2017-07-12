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
validateattributes(valOnQuad, {'cell'}, {'size', [2 1]}, mfilename, 'valOnQuad');
validateattributes(valOnQuad{1}, {'numeric'}, {'size', [R K 3]}, mfilename, 'valOnQuad{1}');
validateattributes(valOnQuad{2}, {'numeric'}, {'size', [R K 3]}, mfilename, 'valOnQuad{2}');

% Assemble matrix
ret = sparse(K*N, KEdge*Nmu);
for n = 1 : 3
  RknTimesVal = sparse(K*N, Nmu);
  markAreaE0T = markE0T(:, n) .* g.areaE0T(:, n);
  for l = 1 : 2
    markAreaSideE0T = markAreaE0T .*  g.markSideE0T(:, n, l);
    for m = 1 : 2
      Rkn = markAreaSideE0T .* g.nuE0T(:, n, m);
      for r = 1 : R
        RknTimesVal = RknTimesVal + kron(Rkn .* valOnQuad{m}(r, :, n)', refEdgePhiIntMuOnQuad{l}(:, :, n, r));
%         ret = ret + ...
%             kron( sparse( 1:K, g.E0T(:, n), ones(K, 1) .* valOnQuad{m}( r, :, n)' .* Rkn , K, KEdge ), ...
%             refEdgePhiIntMuOnQuad{l}( :, :, n, r) );
      end % for r
    end % for m
  end % for l
  ret = ret + kronVec(sparse(1:K, g.E0T(:, n), ones(K, 1), K, KEdge), RknTimesVal);
end % for n
end % function
