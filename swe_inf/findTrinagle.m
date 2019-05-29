function triangles = findTrinagle(g, x1, x2)
%% ATTENTION: NOT OPTIMIZED/VECTORIZED. 
%% This method is intended to be used at the beginnig of a simulation only. 
%% seeks triangles in which a certain pair of coordinates is located via barycentric coordinates
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