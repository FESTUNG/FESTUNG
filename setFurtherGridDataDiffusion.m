function g = setFurtherGridDataDiffusion(g, markE0Tint)
g.areaE0TnuE0T = cell(3,2);
g.RknmarkE0TintnuE0T = cell(3,2);
g.markE0TE0TtimesRkn = cell(3,3,2);
g.markE0TE0TQknnuE0T = cell(3,3,2);
for nn = 1 : 3
  Rkn = 0.5 * g.areaE0T(:,nn);
  for np = 1 : 3
    for m = 1 : 2
      g.markE0TE0TtimesRkn{nn,np,m} = bsxfun(@times, g.markE0TE0T{nn,np}, Rkn.*g.nuE0T(:,nn,m));
      g.markE0TE0TQknnuE0T{nn,np,m} = bsxfun(@times, g.markE0TE0T{nn,np}, Rkn.*g.nuE0T(:,nn,m));
    end % for
  end % for
  for m = 1 : 2
    g.areaE0TnuE0T{nn,m} = g.areaE0T(:,nn).*g.nuE0T(:,nn,m);
    g.RknmarkE0TintnuE0T{nn,m} = Rkn.*markE0Tint(:, nn).*g.nuE0T(:,nn,m);
  end % for
end % for
end % function