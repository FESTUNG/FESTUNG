%%TODO 
%
% I precompute matrices G_bar that allows for an 'easy' evaluation of 
% u_{m} phi_{kj} \partial_{x_m} phi_{ki}
%
function ret = integrateRefElemDphiPhiPerQuad(N, basesOnQuad, qOrd)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
if nargin < 3, p = (sqrt(8*N+1)-3)/2;  qOrd = max(2*p, 1);  end
[~, ~, W] = quadRule2D(qOrd); R = length(W);
ret = { zeros(N, N, R); zeros(N, N, R) };
for i = 1 : N
  for j = 1 : N
    for m = 1 : 2
      ret{m}(i, j, :) =  W(:) .* basesOnQuad.phi2D{qOrd}(:, j) .* basesOnQuad.gradPhi2D{qOrd}(:, i, m);
    end % for
  end % for
end % for
end % function
