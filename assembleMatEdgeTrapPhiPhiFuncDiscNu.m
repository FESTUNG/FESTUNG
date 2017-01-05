function ret = assembleMatEdgeTrapPhiPhiFuncDiscNu(g, markE0Tint, refEdgePhiIntPhiIntPhiInt, refEdgePhiIntPhiExtPhiExt, dataDisc)

if iscell(dataDisc)
  if isequal(size(dataDisc), [2 2])
    ret = assembleMatEdgeTrapPhiPhiFuncDiscMatrixNu(g, markE0Tint, refEdgePhiIntPhiIntPhiInt, refEdgePhiIntPhiExtPhiExt, dataDisc);
  elseif isvector(dataDisc) && numel(dataDisc) == 2
    ret = assembleMatEdgeTrapPhiPhiFuncDiscVectorNu(g, markE0Tint, refEdgePhiIntPhiIntPhiInt, refEdgePhiIntPhiExtPhiExt, dataDisc);
  else
    error('dataDisc must be a KxN-matrix or a 2x1-cell or a 2x2-cell of such matrices.')
  end % if
else
  ret = assembleMatEdgeTrapPhiPhiFuncDiscScalarNu(g, markE0Tint, refEdgePhiIntPhiIntPhiInt, refEdgePhiIntPhiExtPhiExt, dataDisc);
end % if
end % function

function ret = assembleMatEdgeTrapPhiPhiFuncDiscMatrixNu(g, markE0Tint, refEdgePhiIntPhiIntPhiInt, refEdgePhiIntPhiExtPhiExt, dataDisc)
[K, N] = size(dataDisc{1,1});
ret = { sparse(K*N, K*N), sparse(K*N, K*N) };
for n = 1 : 4
  for r = 1 : 2
    areaNuE0Tint = 0.5 * g.areaE0T(:,n) .* g.nuE0T(:,n,r) .* markE0Tint(:,n);
    for m = 1 : 2
      for l = 1 : N
        ret{m} = ret{m} + kron(spdiags(areaNuE0Tint .* dataDisc{r,m}(:,l), 0, K, K), ...
                               refEdgePhiIntPhiIntPhiInt(:,:,l,n)) + ...
                          kron(g.markE0TE0T{n} .* (areaNuE0Tint * dataDisc{r,m}(:,l)'), ...
                               refEdgePhiIntPhiExtPhiExt(:,:,l,n));
      end % for l
    end  % for m
  end % for r
end  % for n
end  % function

function ret = assembleMatEdgeTrapPhiPhiFuncDiscVectorNu(g, markE0Tint, refEdgePhiIntPhiIntPhiInt, refEdgePhiIntPhiExtPhiExt, dataDisc)
[K, N] = size(dataDisc{1});
ret = { sparse(K*N, K*N), sparse(K*N, K*N) };
for n = 1 : 4
  for m = 1 : 2
    areaNuE0Tint = 0.5 * g.areaE0T(:,n) .* g.nuE0T(:,n,m) .* markE0Tint(:,n);
    for l = 1 : N
      ret{m} = ret{m} + kron(spdiags(areaNuE0Tint .* dataDisc{m}(:,l), 0, K, K), ...
                             refEdgePhiIntPhiIntPhiInt(:,:,l,n)) + ...
                        kron(g.markE0TE0T{n} .* (areaNuE0Tint * dataDisc{m}(:,l)'), ...
                             refEdgePhiIntPhiExtPhiExt(:,:,l,n));
    end % for l
  end % for m
end  % for n
end  % function

function ret = assembleMatEdgeTrapPhiPhiFuncDiscScalarNu(g, markE0Tint, refEdgePhiIntPhiIntPhiInt, refEdgePhiIntPhiExtPhiExt, dataDisc)
[K, N] = size(dataDisc);
ret = { sparse(K*N, K*N), sparse(K*N, K*N) };
for n = 1 : 4
  for m = 1 : 2
    areaNuE0Tint = 0.5 * g.areaE0T(:,n) .* g.nuE0T(:,n,m) .* markE0Tint(:,n);
    for l = 1 : N
      ret{m} = ret{m} + kron(spdiags(areaNuE0Tint .* dataDisc(:,l), 0, K, K), ...
                             refEdgePhiIntPhiIntPhiInt(:,:,l,n)) + ...
                        kron(g.markE0TE0T{n} .* (areaNuE0Tint * dataDisc(:,l)'), ...
                             refEdgePhiIntPhiExtPhiExt(:,:,l,n));
    end % for l
  end % for m
end  % for n
end  % function