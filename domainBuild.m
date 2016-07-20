function g = domainBuild(coordV, V0T, OSNodes, LNodes)
%% builds a grid specified by the vertex coordinates and the grid connectivity
g  = generateGridData(coordV, V0T);
g.idE = zeros(g.numE, 1);
sub = 1; % trick
for i = 1 : length(LNodes)
  ind = g.V2E( sub2ind( size(g.V2E), LNodes{i}, [LNodes{i}(2:end); LNodes{i}(1)] ) );
  if i == 12
    sub = 0;
	end
	ind(1:end-sub)
  g.idE(ind(1:end-sub)) = i;
end % for
for i = 1 : length(OSNodes)
  ind = g.V2E( sub2ind( size(g.V2E), OSNodes{i}, [OSNodes{i}(2:end); OSNodes{i}(1)] ) );
  if i == 1
    sub = 1;
  end
  g.idE( ind(1:end-sub) ) = i + length(LNodes);
end % for
%% end TODO
g.idE0T = g.idE(g.E0T); % local edge IDs
end % function