function g = domainHierarchy(X1, X2, initialSize, refinement)
gd = [2; length(X1(:)); X1(:); X2(:)]; % geometry description
sf = 'polygon';                        % set formula
ns = double('polygon')';               % name space
aux = decsg(gd,sf,ns);
[p, e, t] = initmesh(aux, 'Hmax', initialSize);
for i = 1:refinement
  [p, e, t] = refinemesh(aux,p,e,t);
end % for
  g  = generateGridData(p', t(1:3, :)');
  g.idE = zeros(g.numE, 1);
  g.idE(g.V2E(sub2ind(size(g.V2E),e(1,:),e(2,:)))) = e(5,:);
  g.idE0T = g.idE(g.E0T); % local edge IDs
end % function
