function g = setGridData(g, NEDGES, NEDNO, NEDEL, NELED, NLEDN, NRAEDN, NRIEDN, NSEDN)

%% TODO nuE0T


% sollte gehen
g.V2T  = sparse(g.V0T(:, [1 2 3 1 2 3 1 2 3]), g.V0T(:, [2 3 1 2 3 1 2 3 1]), ...
  [(1:g.numT)',zeros(g.numT,3),(1:g.numT)',zeros(g.numT,3),(1:g.numT)'],g.numV,g.numV);

% klappt das???
% the following implicitely defines the edge numbers
[r, c] = find(triu(g.V2T + g.V2T'));
g.V2E = sparse(r, c, 1 : size(r, 1), g.numV, g.numV);
g.V2E = g.V2E + g.V2E';

g.V0E  = NEDNO ;
g.T0E  = NEDEL ;
g.numE = NEDGES;

vecE = g.coordV(g.V0E(:, 2), :) - g.coordV(g.V0E(:, 1), :);
g.areaE = (vecE(:, 1).^2 + vecE(:, 2).^2).^(1/2);

g.nuE = vecE * [0,-1; 1,0] ./ g.areaE(:, [1, 1]);

g.areaT = ...
  ( g.coordV(g.V0T(:,1),1).*g.coordV(g.V0T(:,2),2) + g.coordV(g.V0T(:,2),1).*g.coordV(g.V0T(:,3),2) ...
  + g.coordV(g.V0T(:,3),1).*g.coordV(g.V0T(:,1),2) - g.coordV(g.V0T(:,1),1).*g.coordV(g.V0T(:,3),2) ...
  - g.coordV(g.V0T(:,2),1).*g.coordV(g.V0T(:,1),2) - g.coordV(g.V0T(:,3),1).*g.coordV(g.V0T(:,2),2) )/2;
g.baryT = (g.coordV(g.V0T(:,1),:)+g.coordV(g.V0T(:,2),:)+g.coordV(g.V0T(:,3),:))/3;
g.E0T = NELED;
g.areaE0T = g.areaE(g.E0T);
sigE0T = 1-2*(bsxfun(@eq, reshape(g.T0E(g.E0T,2),g.numT,3),(1:g.numT)'));

g.baryE = 0.5 * (g.coordV(g.V0E(:, 1), :) + g.coordV(g.V0E(:, 2), :));

for n = 1 : 3
  g.coordV0T(:, n, 1) = g.coordV(g.V0T(:, n), 1)';  g.coordV0T(:, n, 2) = g.coordV(g.V0T(:, n), 2)';
  g.baryE0T(:, n, 1) = g.baryE(g.E0T(:, n), 1)';    g.baryE0T(:, n, 2) = g.baryE(g.E0T(:, n), 2)';
  
	%% klappt das?
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
%% klappt das?
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
g.idE = zeros(g.numE, 1);
g.idE(NLEDN ) = 1;
g.idE(NRAEDN) = 2;
g.idE(NRIEDN) = 3;
g.idE(NSEDN ) = 4;
g.idE0T = g.idE(g.E0T); % local edge IDs
end % function