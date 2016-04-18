function triangles = findTrinagle(g, x1, x2)
%% ATTENTION: NOT OPTIMIZED/VECTORIZED. 
%% This method is intended to be used at the beginnig of a simulation only. 
%% uses barycentric coordinates in order to seek triangles in which a certain pair of coordinates is located
%% since these coordinates might be located on an edge or vertex every triangle has to be checked. 
%% The method returns all triangles that are not disjoint with a small ball about the coordinat pair.
K = g.numT;
lambda = zeros(K, 3);
triangles = [];
for k = 1 : K
  lambda(k,:) = [g.coordV0T(k,:,1); g.coordV0T(k,:,2); 1 1 1] \ [x1; x2; 1];
  if lambda(k,:) >= -10^-15*ones(1,3);
    triangles = [triangles, k];
  end % if
end % for
end % function
