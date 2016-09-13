function ret = assembleMatElemTrapDphiPhiFuncDisc(g , refElemDphiPhiPhi, dataDisc)
[K, N] = size(dataDisc{1,1});
ret = { sparse(K*N, K*N), sparse(K*N, K*N) };
for m = 1 : 2
  for r = 1 : 2
    for l = 1 : N
      for s = 1 : 3
        ret{m} = ret{m} + ...
                   kron(spdiags(dataDisc{r,m}(:,l) .* g.J0T{s}(:,3-r,3-r), 0, K, K), refElemDphiPhiPhi{s}(:,:,l,  r)) - ...
                   kron(spdiags(dataDisc{r,m}(:,l) .* g.J0T{s}(:,3-r,  r), 0, K, K), refElemDphiPhiPhi{s}(:,:,l,3-r));
      end % for s
    end % for l
  end % for r
end % for m
end  % function