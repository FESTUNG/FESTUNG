function ret = assembleMatElemTrapDphiPhiFuncDisc(g, refElemDphiPhiPhi, dataDisc)
if iscell(dataDisc)
  if isequal(size(dataDisc), [2 2])
    ret = assembleMatElemTrapDphiPhiFuncDiscMatrix(g, refElemDphiPhiPhi, dataDisc);
  elseif isvector(dataDisc) && numel(dataDisc) == 2
    ret = assembleMatElemTrapDphiPhiFuncDiscVector(g, refElemDphiPhiPhi, dataDisc);
  else
    error('dataDisc must be a KxN-matrix or a 2x1-cell or a 2x2-cell of such matrices.')
  end % if
else
  ret = assembleMatElemTrapDphiPhiFuncDiscScalar(g, refElemDphiPhiPhi, dataDisc);
end % if
end % function

function ret = assembleMatElemTrapDphiPhiFuncDiscMatrix(g, refElemDphiPhiPhi, dataDisc)
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

function ret = assembleMatElemTrapDphiPhiFuncDiscVector(g, refElemDphiPhiPhi, dataDisc)
[K, N] = size(dataDisc{1});
ret = { sparse(K*N, K*N), sparse(K*N, K*N) };
for m = 1 : 2
  for l = 1 : N
    for s = 1 : 3
      ret{m} = ret{m} + ...
                 kron(spdiags(dataDisc{m}(:,l) .* g.J0T{s}(:,3-m,3-m), 0, K, K), refElemDphiPhiPhi{s}(:,:,l,  m)) - ...
                 kron(spdiags(dataDisc{m}(:,l) .* g.J0T{s}(:,3-m,  m), 0, K, K), refElemDphiPhiPhi{s}(:,:,l,3-m));
    end % for s
  end % for l
end % for m
end  % function

function ret = assembleMatElemTrapDphiPhiFuncDiscScalar(g, refElemDphiPhiPhi, dataDisc)
[K, N] = size(dataDisc);
ret = { sparse(K*N, K*N), sparse(K*N, K*N) };
for m = 1 : 2
  for l = 1 : N
    for s = 1 : 3
      ret{m} = ret{m} + ...
                 kron(spdiags(dataDisc(:,l) .* g.J0T{s}(:,3-m,3-m), 0, K, K), refElemDphiPhiPhi{s}(:,:,l,  m)) - ...
                 kron(spdiags(dataDisc(:,l) .* g.J0T{s}(:,3-m,  m), 0, K, K), refElemDphiPhiPhi{s}(:,:,l,3-m));
    end % for s
  end % for l
end % for m
end  % function