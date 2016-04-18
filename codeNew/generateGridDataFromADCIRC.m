function [g, NLEDN, NRAEDN, NRIEDN, NSEDN] = generateGridDataFromADCIRC(g, NEDGES, NEDNO, NLEDN, NRAEDN, NRIEDN, NSEDN)

g.V2T  = sparse(g.V0T(:, [1 2 3 1 2 3 1 2 3]), g.V0T(:, [2 3 1 2 3 1 2 3 1]), ...
  [(1:g.numT)',zeros(g.numT,3),(1:g.numT)',zeros(g.numT,3),(1:g.numT)'],g.numV,g.numV);
% the following implicitely defines the edge numbers
[r, c] = find(triu(g.V2T + g.V2T'));
g.V2E = sparse(r, c, 1 : size(r, 1), g.numV, g.numV);
g.V2E = g.V2E + g.V2E';

idxE = full(g.V2E(sub2ind([g.numV,g.numV],g.V0T(end:-1:1,[1,2,3]),g.V0T(end:-1:1,[2,3,1]))))';
g.V0E(idxE(:), 1) = reshape(g.V0T(end:-1:1, [1,2,3])', 3*g.numT, 1);
g.V0E(idxE(:), 2) = reshape(g.V0T(end:-1:1, [2,3,1])', 3*g.numT, 1);

g.T0E(idxE(:), 1) = reshape(full(g.V2T(sub2ind([g.numV,g.numV], ...
  g.V0T(end:-1:1,[1,2,3]), g.V0T(end:-1:1,[2,3,1]))))', 3*g.numT, 1);
g.T0E(idxE(:), 2) = reshape(full(g.V2T(sub2ind([g.numV,g.numV], ...
  g.V0T(end:-1:1,[2,3,1]), g.V0T(end:-1:1,[1,2,3]))))', 3*g.numT, 1);

g.numE = NEDGES;

vecE = g.coordV(g.V0E(:, 2), :) - g.coordV(g.V0E(:, 1), :);
g.areaE = (vecE(:, 1).^2 + vecE(:, 2).^2).^(1/2);

g.nuE = vecE * [0,-1; 1,0] ./ g.areaE(:, [1, 1]);

g.areaT = ...
  ( g.coordV(g.V0T(:,1),1).*g.coordV(g.V0T(:,2),2) + g.coordV(g.V0T(:,2),1).*g.coordV(g.V0T(:,3),2) ...
  + g.coordV(g.V0T(:,3),1).*g.coordV(g.V0T(:,1),2) - g.coordV(g.V0T(:,1),1).*g.coordV(g.V0T(:,3),2) ...
  - g.coordV(g.V0T(:,2),1).*g.coordV(g.V0T(:,1),2) - g.coordV(g.V0T(:,3),1).*g.coordV(g.V0T(:,2),2) )/2;
g.baryT = (g.coordV(g.V0T(:,1),:)+g.coordV(g.V0T(:,2),:)+g.coordV(g.V0T(:,3),:))/3;

g.E0T = full(g.V2E(sub2ind([g.numV,g.numV],g.V0T(:,[2,3,1]),g.V0T(:,[3,1,2]))));

g.areaE0T = g.areaE(g.E0T);
sigE0T = 1-2*(bsxfun(@eq, reshape(g.T0E(g.E0T,2),g.numT,3),(1:g.numT)'));

g.baryE = 0.5 * (g.coordV(g.V0E(:, 1), :) + g.coordV(g.V0E(:, 2), :));

for n = 1 : 3
  g.coordV0T(:, n, 1) = g.coordV(g.V0T(:, n), 1)';  g.coordV0T(:, n, 2) = g.coordV(g.V0T(:, n), 2)';
  g.baryE0T(:, n, 1) = g.baryE(g.E0T(:, n), 1)';    g.baryE0T(:, n, 2) = g.baryE(g.E0T(:, n), 2)';
  
	g.nuE0T(:, n, 1) = g.nuE(g.E0T(:, n), 1)';        g.nuE0T(:, n, 2) = g.nuE(g.E0T(:, n), 2)';
  Tn = sigE0T(:, n) == 1;  Tp = ~Tn;
  g.E0E(g.E0T(Tn, n), 1) = n;                       g.E0E(g.E0T(Tp, n), 2) = n;
	g.nuE0T(:, n, 1) = g.nuE(g.E0T(:, n), 1)'.* sigE0T(:, n)';
	g.nuE0T(:, n, 2) = g.nuE(g.E0T(:, n), 2)'.* sigE0T(:, n)';
end % for
for m = 1 : 2
  g.B(:, m, 1) = g.coordV0T(:, 2, m) - g.coordV0T(:, 1, m);
  g.B(:, m, 2) = g.coordV0T(:, 3, m) - g.coordV0T(:, 1, m);
end % for
markEint = g.E0E(:, 2) ~= 0; % mark interior edges
g.markE0TE0T = cell(3, 3);
for nn = 1 : 3
  for np = 1 : 3
    g.markE0TE0T{nn,np} = sparse(g.numT, g.numT);
    markEn = g.E0E(:, 1) == nn;  markEp = g.E0E(:, 2) == np;
    idx = markEn & markEp & markEint;
    g.markE0TE0T{nn, np}(sub2ind([g.numT, g.numT], g.T0E(idx, 1), g.T0E(idx, 2))) = 1;
    markEn = g.E0E(:, 2) == nn;  markEp = g.E0E(:, 1) == np;
    idx = markEn & markEp & markEint;
    g.markE0TE0T{nn, np}(sub2ind([g.numT, g.numT], g.T0E(idx, 2), g.T0E(idx, 1))) = 1;
  end % for
end % for

% changed 18.03.16 to make ADCIRC grids consistent to FESTUNG
for k = 1:length(NLEDN)
	[inc, NLEDN(k)] = ismember(NEDNO(NLEDN(k),:), g.V0E, 'rows');
	if ~inc
		error(['The constructed grid has no land edge conecting vertices ' mat2str(g.V0E(NLEDN(k),1)) ' and ' ...
						mat2str(g.V0E(NLEDN(k),1)) ', as in the ADCIRC grid. Check if variables are read in correctly.']);
	end % if
end % for
for k = 1:length(NRAEDN)
	[inc, NRAEDN(k)] = ismember(NEDNO(NRAEDN(k),:), g.V0E, 'rows');
	if ~inc
		error(['The constructed grid has no radiation edge conecting vertices ' mat2str(g.V0E(NRAEDN(k),1)) ' and ' ...
						mat2str(g.V0E(NRAEDN(k),1)) ', as in the ADCIRC grid. Check if variables are read in correctly.']);
	end % if
end % for
for k = 1:length(NRIEDN)
	[inc, NRIEDN(k)] = ismember(NEDNO(NRIEDN(k),:), g.V0E, 'rows');
	if ~inc
		error(['The constructed grid has no river edge conecting vertices ' mat2str(g.V0E(NRIEDN(k),1)) ' and ' ...
						mat2str(g.V0E(NRIEDN(k),1)) ', as in the ADCIRC grid. Check if variables are read in correctly.']);
	end % if
end % for
for k = 1:length(NSEDN)
	[inc, NSEDN(k)] = ismember(NEDNO(NSEDN(k),:), g.V0E, 'rows');
	if ~inc
		error(['The constructed grid has no open sea edge conecting vertices ' mat2str(g.V0E(NSEDN(k),1)) ' and ' ...
						mat2str(g.V0E(NSEDN(k),1)) ', as in the ADCIRC grid. Check if variables are read in correctly.']);
	end % if
end % for

end % function