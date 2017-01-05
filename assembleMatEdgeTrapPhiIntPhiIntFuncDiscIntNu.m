function ret = assembleMatEdgeTrapPhiIntPhiIntFuncDiscIntNu(g, markE0Tbdr, refEdgePhiIntPhiIntPhiInt, dataDisc)
if iscell(dataDisc)
  if isequal(size(dataDisc), [2 2])
    ret = assembleMatEdgeTrapPhiIntPhiIntFuncDiscMatrixIntNu(g, markE0Tbdr, refEdgePhiIntPhiIntPhiInt, dataDisc);
  elseif isvector(dataDisc) && numel(dataDisc) == 2
    ret = assembleMatEdgeTrapPhiIntPhiIntFuncDiscVectorIntNu(g, markE0Tbdr, refEdgePhiIntPhiIntPhiInt, dataDisc);
  else
    error('dataDisc must be a KxN-matrix or a 2x1-cell or a 2x2-cell of such matrices.')
  end % if
else
  ret = assembleMatEdgeTrapPhiIntPhiIntFuncDiscScalarIntNu(g, markE0Tbdr, refEdgePhiIntPhiIntPhiInt, dataDisc);
end % if
end % function

function ret = assembleMatEdgeTrapPhiIntPhiIntFuncDiscMatrixIntNu(g, markE0Tbdr, refEdgePhiIntPhiIntPhiInt, dataDisc)
[K, N] = size(dataDisc{1,1});
ret = { sparse(K*N, K*N), sparse(K*N, K*N) };
for n = 1 : 4
  for r = 1 : 2
    areaNuE0Tbdr = g.areaE0T(:,n) .* g.nuE0T(:,n,r) .* markE0Tbdr(:,n);
    for m = 1 : 2
      for l = 1 : N
        ret{m} = ret{m} + kron(spdiags(areaNuE0Tbdr .* dataDisc{r,m}(:,l), 0, K, K), ...
                               refEdgePhiIntPhiIntPhiInt(:,:,l,n));
      end % for l
    end  % for m
  end % for r
end  % for n
end  % function

function ret = assembleMatEdgeTrapPhiIntPhiIntFuncDiscVectorIntNu(g, markE0Tbdr, refEdgePhiIntPhiIntPhiInt, dataDisc)
[K, N] = size(dataDisc{1});
ret = { sparse(K*N, K*N), sparse(K*N, K*N) };
for n = 1 : 4
  for m = 1 : 2
    areaNuE0Tbdr = g.areaE0T(:,n) .* g.nuE0T(:,n,m) .* markE0Tbdr(:,n);
    for l = 1 : N
      ret{m} = ret{m} + kron(spdiags(areaNuE0Tbdr .* dataDisc{m}(:,l), 0, K, K), ...
                             refEdgePhiIntPhiIntPhiInt(:,:,l,n));
    end % for l
  end  % for m
end  % for n
end  % function

function ret = assembleMatEdgeTrapPhiIntPhiIntFuncDiscScalarIntNu(g, markE0Tbdr, refEdgePhiIntPhiIntPhiInt, dataDisc)
[K, N] = size(dataDisc);
ret = { sparse(K*N, K*N), sparse(K*N, K*N) };
for n = 1 : 4
  for m = 1 : 2
    areaNuE0Tbdr = g.areaE0T(:,n) .* g.nuE0T(:,n,m) .* markE0Tbdr(:,n);
    for l = 1 : N
      ret{m} = ret{m} + kron(spdiags(areaNuE0Tbdr .* dataDisc(:,l), 0, K, K), ...
                             refEdgePhiIntPhiIntPhiInt(:,:,l,n));
    end % for l
  end  % for m
end  % for n
end  % function