function g = refineGrid(g)
K = g.numT;
vertCtr = g.numV;
edges = zeros(g.numE,1);
V0T = zeros(4*K,3);
coordV = [g.coordV; zeros(g.numE,2)];
idE0T = zeros(4*K,3);
bdrTypes = max(g.idE);
for k = 1 : K
	% old vertices
	V0T(4*(k-1)+2,1) = g.V0T(k,1);
	V0T(4*(k-1)+3,2) = g.V0T(k,2);
	V0T(4*(k-1)+4,3) = g.V0T(k,3);
	% local edge 1
	edge = g.E0T(k,1);
	if edges(edge) ~= 0
		V0T(4*(k-1)+1,1) = edges(edge);
		V0T(4*(k-1)+3,3) = edges(edge);
		V0T(4*(k-1)+4,2) = edges(edge);
	else
		vertCtr	= vertCtr+1;
		edges(edge)	= vertCtr;
		coordV(vertCtr,:) = 0.5*(coordV(g.V0E(edge,1),:)+coordV(g.V0E(edge,2),:));
		V0T(4*(k-1)+1,1) = vertCtr;
		V0T(4*(k-1)+3,3) = vertCtr;
		V0T(4*(k-1)+4,2) = vertCtr;
		% boundary edges
		if g.idE(edge) ~= 0
			idE0T(4*(k-1)+3,1) = g.idE(edge);
			idE0T(4*(k-1)+4,1) = g.idE(edge);
		end % if
	end % if
	% local edge 2
	edge = g.E0T(k,2);
	if edges(edge) ~= 0
		V0T(4*(k-1)+1 ,2) = edges(edge);
		V0T(4*(k-1)+2 ,3) = edges(edge);
		V0T(4*(k-1)+4 ,1) = edges(edge);
	else
		vertCtr	= vertCtr+1;
		edges(edge)	= vertCtr;
		coordV(vertCtr,:)	= 0.5*(coordV(g.V0E(edge,1),:)+coordV(g.V0E(edge,2),:));
		V0T(4*(k-1)+1,2) = vertCtr;
		V0T(4*(k-1)+2,3) = vertCtr;
		V0T(4*(k-1)+4,1) = vertCtr;
		% boundary edges
		if g.idE(edge) ~= 0
			idE0T(4*(k-1)+2,2) = g.idE(edge);
			idE0T(4*(k-1)+4,2) = g.idE(edge);
		end % if
	end % if
	% local edge 3
	edge = g.E0T(k,3);
	if edges(edge) ~= 0
		V0T(4*(k-1)+1,3) = edges(edge);
		V0T(4*(k-1)+2,2) = edges(edge);
		V0T(4*(k-1)+3,1) = edges(edge);
	else
		vertCtr	= vertCtr+1;
		edges(edge)	= vertCtr;
		coordV(vertCtr,:)	= 0.5*(coordV(g.V0E(edge,1),:)+coordV(g.V0E(edge,2),:));
		V0T(4*(k-1)+1,3) = vertCtr;
		V0T(4*(k-1)+2,2) = vertCtr;
		V0T(4*(k-1)+3,1) = vertCtr;
		% boundary edges
		if g.idE(edge) ~= 0
			idE0T(4*(k-1)+2,3) = g.idE(edge);
			idE0T(4*(k-1)+3,3) = g.idE(edge);
		end % if
	end % if
end % for
g = generateGridData(coordV,V0T);
g.idE = zeros(g.numE,1);
for b = 1:bdrTypes
  g.idE(g.E0T(idE0T==b)) = b;
end % for
g.idE0T = idE0T;
end % function