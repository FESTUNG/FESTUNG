function globG = assembleMatElemTrapDphiPhiFuncDisc(g , hatG, KDisc)
[K, N] = size(KDisc{1,1});
globG = { sparse(K*N, K*N), sparse(K*N, K*N) };
for m = 1 : 2
  for l = 1 : N
    for r = 1 : 2
      globG{m} = globG{m} + ...
                 kron(spdiags(KDisc{r,m}(:,l) .* g.J0T{s}(:,3-r,3-r), 0, K, K), hatG{s}(:,:,l,  r)) - ...
                 kron(spdiags(KDisc{r,m}(:,l) .* g.J0T{s}(:,3-r,  r), 0, K, K), hatG{s}(:,:,l,3-r));
    end % for r
  end % for l
end % for m
end  % function