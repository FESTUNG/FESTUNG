% TODO
function ret = integrateRefEdgePhiIntMuPerQuad(N, basesOnQuad, qOrd)
validateattributes(N, {'numeric'}, {'numel', 2}, mfilename, 'N')
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
if nargin < 3, p = (sqrt(8*N+1)-3)/2;  qOrd = 2*p+1; end
[~, W] = quadRule1D(qOrd); R = length(W);
ret = { zeros(N(1), N(2), 3, R); zeros(N(1), N(2), 3, R) };
for n = 1 : 3 
  for i = 1 : N(1)
    for j = 1 : N(2)
      ret{1}(i, j, n, :) = W(:) .* basesOnQuad.phi1D{qOrd}(:, i, n) ...
                              .* basesOnQuad.mu{qOrd}(:, j);
      ret{2}(i, j, n, :) = W(:) .* basesOnQuad.phi1D{qOrd}(:, i, n) ...
                              .* basesOnQuad.thetaMu{qOrd}(:, j, 2);
    end % for j
  end % for i
end % for n
end % function
