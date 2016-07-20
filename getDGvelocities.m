function UDG = getDGvelocities(g, cDG, ord, hatM)
%% approximates the velocities by projection of quotients of primary unknowns
global gPhi2D
if length(hatM) == 1
  UDG(:,:,1) = ( cDG(:,:,2) ./ cDG(:,:,1) ) / phi(1,1/3,1/3);
  UDG(:,:,2) = ( cDG(:,:,3) ./ cDG(:,:,1) ) / phi(1,1/3,1/3);
else
  ord = max(ord,1);  [~, ~, W] = quadRule2D(ord);
  K = g.numT; N = size(hatM, 1);
  UDG = zeros(K,N,2);
  UDG(:,:,1) = ( cDG(:,:,2) * gPhi2D{ord}.' ) ./ ( cDG(:,:,1) * gPhi2D{ord}.' ) * (repmat(W.', 1, N) .* gPhi2D{ord}) / hatM;
  UDG(:,:,2) = ( cDG(:,:,3) * gPhi2D{ord}.' ) ./ ( cDG(:,:,1) * gPhi2D{ord}.' ) * (repmat(W.', 1, N) .* gPhi2D{ord}) / hatM;
end % if
end % function
