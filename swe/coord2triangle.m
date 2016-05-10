function ret = coord2triangle(g, coord1, coord2)
numCoord = length(coord1);
% Compute barycentric coordinates for each cartesian coordinate pair
coordBary = zeros(numCoord, g.numT, 3);
coordBary(:,:,1) = bsxfun(@times, (g.coordV0T(:,2,2) - g.coordV0T(:,3,2)).', ...
                          bsxfun(@minus, coord1, g.coordV0T(:,3,1).')) + ...
                   bsxfun(@times, (g.coordV0T(:,3,1) - g.coordV0T(:,2,1)).', ...
                          bsxfun(@minus, coord2, g.coordV0T(:,3,2).'));
coordBary(:,:,2) = bsxfun(@times, (g.coordV0T(:,3,2) - g.coordV0T(:,1,2)).', ...
                          bsxfun(@minus, coord1, g.coordV0T(:,3,1).')) + ...
                   bsxfun(@times, (g.coordV0T(:,1,1) - g.coordV0T(:,3,1)).', ...
                          bsxfun(@minus, coord2, g.coordV0T(:,3,2).'));
coordBary(:,:,1:2) = bsxfun(@times, coordBary(:,:,1:2), .5 ./ g.areaT.');
coordBary(:,:,3) = 1 - coordBary(:,:,1) - coordBary(:,:,2);
% Find triangles with only positive barycentric coordinates
[idxCoord, idxTriangle] = find(coordBary(:,:,1) >= 0 & coordBary(:,:,2) >= 0 & coordBary(:,:,3) >= 0);
% TODO: vectorize this
ret = cell(numCoord,1);
for n = 1 : numCoord
  triangles = idxTriangle(idxCoord == n);
  ret{n} = [triangles, reshape(coordBary(n,triangles,:), length(triangles), 3) ];
end % for
end % function
