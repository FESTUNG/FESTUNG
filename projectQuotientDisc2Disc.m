function ret = projectQuotientDisc2Disc(dataDisc1, dataDisc2, ord, refElemPhiPhi, basesOnQuad)
[~, N] = size(dataDisc1);
validateattributes(dataDisc1, {'numeric'}, {}, mfilename, 'dataDisc1')
validateattributes(dataDisc2, {'numeric'}, {'size', size(dataDisc1)}, mfilename, 'dataDisc2')
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')

if N == 1
  ret = ( dataDisc1 ./ dataDisc2 ) / phi(1,1/3,1/3);
else
  ord = max(ord,1);  [~, ~, W] = quadRule2D(ord);
  ret = ( dataDisc1 * basesOnQuad.phi2D{ord}.' ) ./ ( dataDisc2 * basesOnQuad.phi2D{ord}.' ) * ...
    (repmat(W.', 1, N) .* basesOnQuad.phi2D{ord}) / refElemPhiPhi;
end % if
end % function
