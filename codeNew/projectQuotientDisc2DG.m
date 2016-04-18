function ret = projectQuotientDisc2DG(dataDG1, dataDG2, ord, refElemPhiPhi)
global gPhi2D
[~, N] = size(dataDG1);
if N == 1
  ret = ( dataDG1 ./ dataDG2 ) / phi(1,1/3,1/3);
else
  ord = max(ord,1);  [~, ~, W] = quadRule2D(ord);
  ret = ( dataDG1 * gPhi2D{ord}.' ) ./ ( dataDG2 * gPhi2D{ord}.' ) * (repmat(W.', 1, N) .* gPhi2D{ord}) / refElemPhiPhi;
end % if
end % function
